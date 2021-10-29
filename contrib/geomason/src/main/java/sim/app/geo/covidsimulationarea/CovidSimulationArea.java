/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "LICENSE" for more information
 *
 * $Id$
*/
package sim.app.geo.covidsimulationarea;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.planargraph.Node;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.field.geo.GeomVectorField;
import sim.io.geo.ShapeFileImporter;
import sim.util.geo.GeomPlanarGraph;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 * The simulation core.
 * <p/>
 * The simulation can require a LOT of memory, so make sure the virtual machine
 * has enough. Do this by adding the following to the command line, or by
 * setting up your run configuration in Eclipse to include the VM argument:
 * <p/>
 * -Xmx2048M
 * <p/>
 * With smaller simulations this chunk of memory is obviously not necessary. You
 * can take it down to -Xmx800M or some such. If you get an OutOfMemory error,
 * push it up.
 */
public class CovidSimulationArea extends SimState
{
    private static final long serialVersionUID = 1L;

    int population = 50; // total population initialization

    /**
     * Main function allows simulation to be run in stand-alone, non-GUI mode
     */
    public static void main(String[] args)
    {
        doLoop(CovidSimulationArea.class, args);
        System.exit(0);
    }
    
    public boolean init = true;
    
    public int cycle_count = 0; // used for printing log
    public int total_edges = 0;
    public int total_buildings = 0;
    public int init_infect = 1;
    public GeomVectorField roads = new GeomVectorField();
    public GeomVectorField buildings = new GeomVectorField();
    
    public boolean mitigation_switch = true;
    public int mitigation_threshold = 125;
    public boolean orig_params = true;
    
    // traversable network
    public GeomPlanarGraph network = new GeomPlanarGraph();
    public GeomVectorField junctions = new GeomVectorField();
    
    // mappings for edges
    // mapping between unique edge IDs and edge structures themselves
    public HashMap<Integer, GeomPlanarGraphEdge> idsToEdges =
        new HashMap<Integer, GeomPlanarGraphEdge>();
    public HashMap<GeomPlanarGraphEdge,Integer> edgesToIds =
            new HashMap<GeomPlanarGraphEdge,Integer>();
    public HashMap<Coordinate, GeomPlanarGraphEdge> edgeCoordToEdges = 
        new HashMap<Coordinate, GeomPlanarGraphEdge>();
    
    // mappings for buildings
    public HashMap<Coordinate, GeomPlanarGraphEdge> buildingsToEdges = 
    	new HashMap<Coordinate, GeomPlanarGraphEdge>();      
    public Vector<Coordinate> idsToBuildings = new Vector<Coordinate>();
    public HashMap<Coordinate, Integer> buildingsToVisits = 
            new HashMap<Coordinate, Integer>();
    
    // mappings for areas
    public HashMap<Integer,ArrayList<GeomPlanarGraphEdge>> areasToEdges = 
            new HashMap<Integer,ArrayList<GeomPlanarGraphEdge>>();
    public Vector<Integer> idsToAreas = new Vector<Integer>();
    public HashMap<Coordinate, Integer> buildingsToAreas = 
            new HashMap<Coordinate, Integer>();
    public HashMap<Integer,ArrayList<Coordinate>> areasToBuildings = 
            new HashMap<Integer,ArrayList<Coordinate>>();
    
    // mappings for hotspots
    public HashMap<Coordinate, Coordinate> buildingsToHotspots = 
            new HashMap<Coordinate, Coordinate>();
    public HashMap<Coordinate,ArrayList<Coordinate>> hotspotsToBuildings = 
            new HashMap<Coordinate,ArrayList<Coordinate>>();       
    public HashMap<Coordinate, Double> hotspotsToProbs = 
            new HashMap<Coordinate, Double>();
    public HashMap<Coordinate, Double> hotspotsToRadius = 
            new HashMap<Coordinate, Double>();  
    public HashMap<Integer, ArrayList<Coordinate>> areasToHotspots = 
            new HashMap<Integer, ArrayList<Coordinate>>();
    
    // mappings for types
    public HashMap<Coordinate, String> buildingsToTypes = 
            new HashMap<Coordinate, String>();
    public HashMap<Double,String> cumulativeProbType = new HashMap<Double,String>();
    
    // agents
    public GeomVectorField infected_agents = new GeomVectorField();
    public GeomVectorField healthy_agents = new GeomVectorField();
    ArrayList<Agent> agentList = new ArrayList<Agent>();
    
    // some default parameters
    boolean read_full_map = false; // if true, read the whole map of Seoul.
    double infect_probability = 0.004;  
    public double caution_factor = 1.0;
    
    // variables related to distribution
    // distance
    public ArrayList<Double> xAxisDistance = new ArrayList<>();
    public ArrayList<Double> yAxisDistance = new ArrayList<>();
    
    // mobility
    ArrayList<Double> xAxisLowSpreader = new ArrayList<>();
    ArrayList<Double> yAxisLowSpreader = new ArrayList<>();
    
    ArrayList<Double> xAxisSuperSpreader = new ArrayList<>();
    ArrayList<Double> yAxisSuperSpreader = new ArrayList<>();
    
    // active after symptom onset
    ArrayList<Double> xAxisSymptomMove = new ArrayList<>();
    ArrayList<Double> yAxisSymptomMove = new ArrayList<>();
    
    
    /**
     * Constructor
     */
    public CovidSimulationArea(long seed)
    {
        super(seed);
    }

    /**
     * Initialization
     */
    @Override
    public void start()
    {
        super.start();
        try
        {
        	System.out.println("Beginning log.");
        	// Read population and infect rate
        	Scanner s = new Scanner(new File("data/parameters.txt"));
        	population = s.nextInt();
        	infect_probability = s.nextDouble();
        	caution_factor = s.nextDouble();
        	init_infect = s.nextInt();
        	mitigation_threshold=s.nextInt();
        	s.close();
        	population = 100; // for testing
        	//population=500;
        	System.out.println("rate_infection= " + infect_probability);
        	System.out.println("pop_infection= " + population);
        	System.out.println("caution_factor= " + caution_factor);
        	System.out.println("init_infect= "+init_infect);
            // Read the X and Y axis of the distribution from external files
        	System.out.println("reading distributions...");
            
            Scanner sDistance = new Scanner(new File("data/distr/xAxisDistance.txt"));
            xAxisDistance = new ArrayList<>();
            while(sDistance.hasNext()) {
            	xAxisDistance.add(sDistance.nextDouble());
            }
            sDistance.close();            
            
            Scanner sLowSpreader = new Scanner(new File("data/distr/xAxisLowSpreader.txt"));
            xAxisLowSpreader = new ArrayList<>();
            while(sLowSpreader.hasNext()) {
            	xAxisLowSpreader.add(sLowSpreader.nextDouble());
            }
            sLowSpreader.close();
            
            Scanner sSuperSpreader = new Scanner(new File("data/distr/xAxisSuperSpreader.txt"));
            xAxisSuperSpreader = new ArrayList<>();
            while(sSuperSpreader.hasNext()) {
            	xAxisSuperSpreader.add(sSuperSpreader.nextDouble());
            }
            sSuperSpreader.close();
            
            sDistance = new Scanner(new File("data/distr/yAxisDistance.txt"));
            yAxisDistance = new ArrayList<>();
            while(sDistance.hasNext()) {
            	yAxisDistance.add(sDistance.nextDouble());
            }
            sDistance.close();            
            
            sLowSpreader = new Scanner(new File("data/distr/yAxisLowSpreader.txt"));
            yAxisLowSpreader = new ArrayList<>();
            while(sLowSpreader.hasNext()) {
            	yAxisLowSpreader.add(sLowSpreader.nextDouble());
            }
            sLowSpreader.close();
            
            sSuperSpreader = new Scanner(new File("data/distr/yAxisSuperSpreader.txt"));
            yAxisSuperSpreader = new ArrayList<>();
            while(sSuperSpreader.hasNext()) {
            	yAxisSuperSpreader.add(sSuperSpreader.nextDouble());
            }
            sSuperSpreader.close();
            
            s = new Scanner(new File("data/distr/xAxisSymptomMove.txt"));
            xAxisSymptomMove = new ArrayList<>();
            while(s.hasNext()) {
            	xAxisSymptomMove.add(s.nextDouble());
            }
            s.close();
            
            s = new Scanner(new File("data/distr/yAxisSymptomMove.txt"));
            yAxisSymptomMove = new ArrayList<>();
            while(s.hasNext()) {
            	yAxisSymptomMove.add(s.nextDouble());
            }
            s.close();
            
            Scanner typeProbFile = new Scanner(new File("data/type_probs.txt"));
            Double prob_sum=0.0;
            Double prob;
            String loc_type;
            while(typeProbFile.hasNext()) {
            	loc_type=typeProbFile.next();
            	prob=typeProbFile.nextDouble();
            	prob_sum+=prob;
            	cumulativeProbType.put(prob_sum,loc_type);
            }
            typeProbFile.close();
        	
            // read in the roads to create the transit network
            System.out.println("reading roads layer...");
            URL roadsFile =  Paths.get("data","gangnam_seocho_new","gangnam_seocho_roads.shp").toUri().toURL() ; 
            URL roadsDB =  Paths.get("data","gangnam_seocho_new","gangnam_seocho_roads.dbf").toUri().toURL() ;
            if (read_full_map == true)
            {
            	roadsFile =  Paths.get("data","seoul_roads_split.shp").toUri().toURL() ; 
            	roadsDB =  Paths.get("data","seoul_roads_split.dbf").toUri().toURL() ; 
            }
            ShapeFileImporter.read(roadsFile, roadsDB, roads);
            Envelope MBR = roads.getMBR();

            // read in the tracts to create the background
            System.out.println("reading tracts layer...");
            URL areasFile = Paths.get("data","gangnam_seocho_new","gangnam_seocho_buildings.shp").toUri().toURL() ; 
            URL areasDB = Paths.get("data","gangnam_seocho_new","gangnam_seocho_buildings.dbf").toUri().toURL() ;           
            if (read_full_map == true)
            {
                areasFile = Paths.get("data","seoul_building_area.shp").toUri().toURL() ; 
                areasDB = Paths.get("data","seoul_building_area.dbf").toUri().toURL() ;            	
            }
            ShapeFileImporter.read(areasFile, areasDB, buildings);
            MBR.expandToInclude(buildings.getMBR());

            createNetwork();

            // update so that everyone knows what the standard MBR is
            roads.setMBR(MBR);
            buildings.setMBR(MBR);

            // initialize agents
            populate();
            infected_agents.setMBR(MBR);
            healthy_agents.setMBR(MBR);
            
            // Ensure that the spatial index is updated after all the agents
            // move
            schedule.scheduleRepeating(infected_agents.scheduleSpatialIndexUpdater(), Integer.MAX_VALUE, 1.0);
            schedule.scheduleRepeating(healthy_agents.scheduleSpatialIndexUpdater(), Integer.MAX_VALUE, 1.0);

            /**
             * Steppable that flips Agent paths once everyone reaches their
             * destinations
             */
            Steppable flipper = new Steppable()
            {
                @Override
                public void step(SimState state)
                {
                	// status to sync after each cycle.
                	cycle_count ++;
                	
                	String fileContent=String.valueOf(buildingsToVisits);
                	
                	BufferedWriter fileName;
					try {
						fileName = new BufferedWriter(new FileWriter("data/buildingsToVisits_pop"+String.valueOf(population)+"_rate"+String.valueOf(infect_probability)+"_factor"+caution_factor+".txt"));
						fileName.write(fileContent);
						fileName.close();
					} catch (IOException e) {
						e.printStackTrace();
						System.out.println("error writing to building visit file");
						System.out.println(e);
					}

                	
                	if (cycle_count > 14440) // 50 days.
                	{
                    	System.out.println("<log_fin>cycle="+ cycle_count+", simulated 50 days, exit.");
                    	System.exit(0);
                	}
                	
                	CovidSimulationArea gstate = (CovidSimulationArea) state;

                	// After each cycle, get all agents' location, and compare their destination building.
                	HashMap<Coordinate, ArrayList<Agent>> buildingToAgents =
                	        new HashMap<Coordinate, ArrayList<Agent>>();

                	// Agents to be isolated.
                	//ArrayList<Agent> toBeRemoved = new ArrayList<Agent>();
                	
                    for (Agent a : gstate.agentList)
                    {                    	
                    	if (a.getInBuildingStatus() == true)
                    		// a is inside a building
                    	{
                    		Coordinate curBuilding = a.getDestinationBuilding();
                    		ArrayList<Agent> tmp = buildingToAgents.get(curBuilding);
                    		// check the status of this building: is there any agent also in this building?
                    		if (tmp == null)
                    		{
                    			tmp = new ArrayList<Agent>();
                    		}
                    		tmp.add(a);
                    		buildingToAgents.put(curBuilding, tmp);
                    		
                    		// print log: we find a in the building.
                    		System.out.println("<log> cycle="+ cycle_count + " agent="+ a + " building="+ curBuilding + " infection=" + a.checkInfection());
                    	}
                    }
                    
                    // infection happens in the buildings.
                    Set<Entry<Coordinate, ArrayList<Agent>>> curbuildings = buildingToAgents.entrySet();
                    Iterator<Entry<Coordinate, ArrayList<Agent>>> iterator = curbuildings.iterator();
                    while(iterator.hasNext()) {
                       Map.Entry<Coordinate, ArrayList<Agent>> mentry = (Map.Entry<Coordinate, ArrayList<Agent>>)iterator.next();
                       ArrayList<Agent> curAgents = (ArrayList<Agent>) mentry.getValue();
                       infect(curAgents);                       
                    }
                    
                    // count total number of infection after this cycle.
                    int infectCount = 0;
                    for (Agent a : gstate.agentList)
                    {
                    	if (a.checkInfection() == true)
                    	{
                    		infectCount ++;
                    	}
                    }
                    System.out.println("<log> cycle="+ cycle_count + " infection="+ infectCount);
                    
                    if (mitigation_switch && (infectCount>=mitigation_threshold) && orig_params) {
                    	//change params if imposing stay-at-home order at a certain threshold
                    	Scanner new_params;
						try {
							new_params = new Scanner(new File("data/parameters_new.txt"));
							caution_factor = new_params.nextDouble();
	                    	new_params.close();
						} catch (FileNotFoundException e) {
							e.printStackTrace();
							System.out.println("New parameter file not found");
						}
                    	
                    	
                    	orig_params=false;
                    }
                    
                    // if everyone is infected, then stop simulation.
                    if (infectCount == population)
                    {
                    	System.out.println("<log_fin>Everyone is infected. Simulation finishes!");
                    	System.exit(0);
                    }
                }
            };
            schedule.scheduleRepeating(flipper, 1);

        } catch (FileNotFoundException ex)
        {
            System.out.println("Error: missing required data file");
        } catch (IOException ex)
        {
            Logger.getLogger(CovidSimulationArea.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex)
        {
            Logger.getLogger(CovidSimulationArea.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void infect(ArrayList<Agent> curAgents)
    {
    	int infected = 0;
    	for (Agent a : curAgents)
    	{
    		if (a.checkInfection() == true)
    		{
    			infected ++;
    		}
    	}
    	
    	double pr_infection = 1; 
    	for (int i = 0; i < infected; i++)
    	{
    		pr_infection = (1-infect_probability) * pr_infection;
    	}
    	
    	pr_infection = 1 - pr_infection;
    	
    	if (infected > 0)
    	{
        	for (Agent a : curAgents)
        	{	
        		if (a.checkInfection() == false)
        		{
        			double tmp = ThreadLocalRandom.current().nextDouble();
        			if (tmp < pr_infection)
        			{
	        			a.setInfection();
	        			healthy_agents.removeGeometry(a.getGeometry());
	        			infected_agents.addGeometry(a.getGeometry());
	            		System.out.println("<log-infection> cycle="+ cycle_count + " agent="+ a);
        			}
        		}
        	}
    	}
    }

    /**
     * Create the road network the agents will traverse
     * <p/>
     */
    private void createNetwork()
    {
    	init = true;
        System.out.println("creating network...");
        
        int a_count=0;
    	ArrayList<Double> origin_areas = new ArrayList<Double>();
        Scanner origin_probs_file;
		try {
			origin_probs_file = new Scanner(new File("data/distr/origin_probs.txt")); //just to get number of areas...will replace
			while(origin_probs_file.hasNext())
	    	{
	    		origin_areas.add(Double.valueOf(a_count));
	    		Double orig_prob=origin_probs_file.nextDouble();
	    		a_count++;
	    	}
	    	origin_probs_file.close();
		} catch (FileNotFoundException e1) {
			System.out.println("origin area probability file not found");
			e1.printStackTrace();
			System.exit(1);
		}     	
		
    	int num_areas=a_count;
    	
    	for (int i=0;i<num_areas;i++) {
    		areasToBuildings.put(i,new ArrayList<Coordinate>());
    		areasToEdges.put(i,new ArrayList<GeomPlanarGraphEdge>());
    		areasToHotspots.put(i, new ArrayList<Coordinate>());
    	}
    	
		try {
			Scanner fileHotspotArea;
			fileHotspotArea = new Scanner(new File("data/hotspot_area.txt"));
		
			while(fileHotspotArea.hasNext()) {
				areasToHotspots.get(fileHotspotArea.nextInt()).add(new Coordinate(fileHotspotArea.nextDouble(),fileHotspotArea.nextDouble()));
			}
			fileHotspotArea.close();
        
			Scanner fileHotspotProbs = new Scanner(new File("data/hotspot_prob.txt"));
			while(fileHotspotProbs.hasNext()) {
				hotspotsToProbs.put(new Coordinate(fileHotspotProbs.nextDouble(),fileHotspotProbs.nextDouble()),fileHotspotProbs.nextDouble());
			}
			fileHotspotProbs.close();
        
			Scanner fileHotspotRadius = new Scanner(new File("data/hotspot_radius.txt"));
			while(fileHotspotRadius.hasNext()) {
				Coordinate c = new Coordinate(fileHotspotRadius.nextDouble(),fileHotspotRadius.nextDouble());
				hotspotsToRadius.put(c,fileHotspotRadius.nextDouble());
				hotspotsToBuildings.put(c,new ArrayList<Coordinate>());
			}
			fileHotspotRadius.close();
		} catch (FileNotFoundException e1) {
			System.out.println("no hotspot files");
		}

        network.createFromGeomField(roads);
        int countedge = 0;
        for (Object o : network.getEdges())
        {
            GeomPlanarGraphEdge e = (GeomPlanarGraphEdge) o;
            idsToEdges.put(countedge, e);
            edgesToIds.put(e, countedge);
            countedge ++;

            e.setData(new ArrayList<Agent>());

        	Node curNode = e.getDirEdge(0).getToNode();
            Coordinate curCoord = curNode.getCoordinate();
            
            edgeCoordToEdges.put(curCoord, e);
            curNode = e.getDirEdge(1).getToNode();
            curCoord = curNode.getCoordinate();
            edgeCoordToEdges.put(curCoord, e);
        }

        addIntersectionNodes(network.nodeIterator(), junctions);
        total_edges = idsToEdges.size();

        // dealing with buildings
        
        // create a hashmap, from building's Coordinate to edges(roads).
        for (Object o : buildings.getGeometries())
        {        	
            MasonGeometry b = (MasonGeometry) o;
            String type = b.getStringAttribute("type");
            Coordinate cur_b = b.getGeometry().getCoordinate();
          
            int area_code=b.getIntegerAttribute("area_code");
            
            
            GeomPlanarGraphEdge cur_edge = new GeomPlanarGraphEdge(null);
            double min_dist = 1000000000;
            
            // go through building coordinates to find the closest
            Set<Entry<Coordinate, GeomPlanarGraphEdge>> coords = edgeCoordToEdges.entrySet();
            Iterator<Entry<Coordinate, GeomPlanarGraphEdge>> iterator = coords.iterator();
            while(iterator.hasNext()) {
               Map.Entry<Coordinate, GeomPlanarGraphEdge> mentry = (Map.Entry<Coordinate, GeomPlanarGraphEdge>)iterator.next();
               
               double cur_dist = getCoordDistance(cur_b, mentry.getKey());
               if (cur_dist < min_dist)
               {
            	   min_dist = cur_dist;
            	   cur_edge = mentry.getValue();
               }
            }
            
            
            if (min_dist >= 1000000000)
            {
            	System.exit(0);
            }
            buildingsToVisits.put(cur_b, 0);
            buildingsToAreas.put(cur_b,area_code);
            buildingsToEdges.put(cur_b, cur_edge);
            buildingsToTypes.put(cur_b, type);
            idsToBuildings.add(cur_b);
            areasToBuildings.get(area_code).add(cur_b);
            areasToEdges.get(area_code).add(cur_edge);
            idsToAreas.add(area_code);
            
            ArrayList<Coordinate> hotspots_area = areasToHotspots.get(area_code);
            for (int h=0;h<hotspots_area.size();h++) {
            	Coordinate spot=hotspots_area.get(h);
            	if (getCoordDistance(spot,cur_b) <= hotspotsToRadius.get(spot)) {
            		
            		buildingsToHotspots.put(cur_b, spot);
            		hotspotsToBuildings.get(spot).add(cur_b);
            		
            	}
            }
            
        }

        total_buildings = buildingsToEdges.size();
        System.out.println("create network DONE in total " +buildingsToEdges.size() +" buildings");
    }
    
    /**
     * create an appropriate pop
     * <p/>     
     */
    public void populate()
    {
        try
        {
        	Scanner origin_probs_file= new Scanner(new File("data/distr/origin_probs.txt"));
        	ArrayList<Double> origin_probs = new ArrayList<Double>();
        	int a_count=0;
        	ArrayList<Double> origin_areas = new ArrayList<Double>();
        	while(origin_probs_file.hasNext())
        	{
        		origin_probs.add(origin_probs_file.nextDouble());
        		origin_areas.add(Double.valueOf(a_count));
        		a_count++;
        	}
        	origin_probs_file.close();
        	Scanner origin_superspreader = new Scanner(new File("data/area_superspreader.txt"));
        	ArrayList<Double> superspreader_probs = new ArrayList<Double>();
        	while(origin_superspreader.hasNext())
        	{
        		superspreader_probs.add(origin_superspreader.nextDouble());
        	}
        	int num_areas=a_count;
        	Scanner area_probs = new Scanner(new File("data/distr/area_probs.txt"));
        	ArrayList<ArrayList<Double>> area_probs_all = new ArrayList<ArrayList<Double>>();
        	
        	for (int j=0;j<num_areas;j++)
        	{
        		area_probs_all.add(new ArrayList<Double>());
        		for (int k=0;k<num_areas;k++)
        		{
        			area_probs_all.get(j).add(area_probs.nextDouble());
        		}
        	}
        	
        	area_probs.close();
        	
        	Scanner origin_stay_file = new Scanner(new File("data/distr/area_stay_probs.txt"));
        	ArrayList<Double> origin_stay_probs = new ArrayList<Double>();
        	for (int m=0;m<num_areas;m++)
        	{
        		origin_stay_probs.add(origin_stay_file.nextDouble());
        	}
        	origin_stay_file.close();
        	String fileContent="";
        	for (int i = 0; i < population; i++) 
            {
            	System.out.println("in population "+i);
                
                int startPos = ThreadLocalRandom.current().nextInt(total_edges);
                GeomPlanarGraphEdge startingEdge =
                    idsToEdges.get(startPos);  
                Agent a = new Agent(this, startingEdge);
                a.init = true;
                
                a.setArea(origin_areas, origin_probs,  origin_stay_probs, area_probs_all);
                
                // super spreaders and low spreaders + mobility
                if (ThreadLocalRandom.current().nextDouble() < superspreader_probs.get(a.area) ) 
                {
                	a.superspreader = true;
                	a.setMobility(xAxisSuperSpreader, yAxisSuperSpreader);
                	System.out.println("Agent="+a+" set superspreader");
                }
                else
                {
                	a.setMobility(xAxisLowSpreader, yAxisLowSpreader);
                	System.out.println("Agent="+a+" set lowspreader");
                }
                
                // speed
            	// all these small double numbers are degree related to the earth radius. its unit is not mile or meter.
                if (i % 2 == 0)
                {
                	// walk
                	a.speed = 5*0.0035972;
                	a.move_speed = 5*0.0035972;
                }
                else
                {
                	// car
                	double maxspeed = 0.029977;
                	double minspeed = 0.01124125;
                	double tmp = ThreadLocalRandom.current().nextDouble();
                	a.speed = 5*(minspeed + (maxspeed-minspeed)*tmp);
                	a.move_speed = a.speed;
                }
                
                
                ArrayList<Coordinate> blist = areasToBuildings.get(a.checkArea());
                int b_pos=ThreadLocalRandom.current().nextInt(blist.size());
                startingEdge = buildingsToEdges.get(blist.get(b_pos));
                a.startNode = startingEdge.getDirEdge(0).getFromNode();
                a.homeBldg=blist.get(b_pos);
                GeomPlanarGraphEdge goalEdge = a.chooseDestination();
                a.endNode = goalEdge.getDirEdge(0).getToNode();
                Coordinate startCoord = a.startNode.getCoordinate(); 
                a.updatePosition(startCoord);
                while(!a.start(this)) {
                	b_pos=ThreadLocalRandom.current().nextInt(blist.size());
                    startingEdge = buildingsToEdges.get(blist.get(b_pos));
                      
                    a.startNode = startingEdge.getDirEdge(0).getFromNode(); 
                    a.homeBldg=blist.get(b_pos);
                    goalEdge = a.chooseDestination();
                    a.endNode = goalEdge.getDirEdge(0).getToNode();
                    
                    startCoord = null;
                    startCoord = a.startNode.getCoordinate();
                    a.updatePosition(startCoord);
                }

                MasonGeometry newGeometry = a.getGeometry();
                newGeometry.isMovable = true;
                
                // set infection
//                if (i < init_infection_rate * population)
//                if (i < 1)
                if (i < init_infect)
                {
                	a.setInfection();
                	infected_agents.addGeometry(newGeometry);
                }
                else
                {
                	healthy_agents.addGeometry(newGeometry);
                }
 
                agentList.add(a);
                schedule.scheduleRepeating(a);
                
                init = false;
                a.init = false;
                fileContent=fileContent+String.valueOf(a.getGeometry())+"\t";
                fileContent=fileContent+String.valueOf(a.superspreader)+"\t";
                fileContent=fileContent+String.valueOf(a.getMobility())+"\t";
                fileContent=fileContent+String.valueOf(a.checkInfection())+"\t";
                fileContent=fileContent+String.valueOf(a.checkArea())+"\t";
                fileContent=fileContent+String.valueOf(a.checkAreaStay())+"\t";
                fileContent=fileContent+String.valueOf(a.move_speed)+"\n";
                
                //Path fileName = Path.of("data/area2_profiling_pop5000.txt");                
                //Files.writeString(fileName, fileContent);
                
             }
        	
        } catch (Exception e)
        {
        	System.out.println(e);
            System.out.println("error in populate()");
            System.exit(1);
        }
    }

    /**
     * adds nodes corresponding to road intersections to GeomVectorField
     * <p/>
     * @param nodeIterator  Points to first node
     * @param intersections GeomVectorField containing intersection geometry
     * <p/>
     * Nodes will belong to a planar graph populated from LineString network.
     */
    private void addIntersectionNodes(Iterator<?> nodeIterator,
                                      GeomVectorField intersections)
    {
        GeometryFactory fact = new GeometryFactory();
        Coordinate coord = null;
        Point point = null;

        while (nodeIterator.hasNext())
        {
            Node node = (Node) nodeIterator.next();
            coord = node.getCoordinate();
            point = fact.createPoint(coord);

            junctions.addGeometry(new MasonGeometry(point));
        }
    }
    
    public double getCoordDistance(Coordinate a, Coordinate b) //returns distance in meters
    {
    	double lng1 = a.x;    // longitude
    	double lng2 = b.x;
    	double lat1 = a.y;
    	double lat2 = b.y;
    
    	double x = (lng2 - lng1) * Math.PI * 6.371229 * 1e6 * Math.cos(((lat1 + lat2) / 2) * Math.PI / 180) / 180;
        double y = (lat1 - lat2) * Math.PI * 6.371229 * 1e6 / 180;
   
        return Math.hypot(x, y);
    }

}

