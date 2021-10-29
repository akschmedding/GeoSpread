/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "LICENSE" for more information
 *
 * $Id$
*/
package sim.app.geo.covidsimulationsymptom;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Envelope;
import com.vividsolutions.jts.geom.GeometryFactory;
import com.vividsolutions.jts.geom.Point;
import com.vividsolutions.jts.planargraph.Node;

import java.io.File;
import java.io.FileNotFoundException;
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
public class CovidSimulationSymptom extends SimState
{
    private static final long serialVersionUID = 1L;

    int population = 50; // total population 
    int total_isolated = 0; // agents that are isolated
    /**
     * Main function allows simulation to be run in stand-alone, non-GUI mode
     */
    public static void main(String[] args)
    {
        doLoop(CovidSimulationSymptom.class, args);
        System.exit(0);
    }
    public boolean init = true;
    
    public int cycle_count = 0; // used for printing log
    public int total_edges = 0;
    public int total_buildings = 0;
    double init_infection_rate = 0.01; // Percentage of infected agents at the beginning. Currently this is not used. We start from 1 infected agent.
    
    public GeomVectorField roads = new GeomVectorField();
    public GeomVectorField buildings = new GeomVectorField();
    
    // traversable network
    public GeomPlanarGraph network = new GeomPlanarGraph();
    public GeomVectorField junctions = new GeomVectorField();
    
    // mapping between unique edge IDs and edge structures themselves
    public HashMap<Integer, GeomPlanarGraphEdge> idsToEdges =
        new HashMap<Integer, GeomPlanarGraphEdge>();
    public HashMap<Coordinate, GeomPlanarGraphEdge> buildingsToEdges = 
    	new HashMap<Coordinate, GeomPlanarGraphEdge>();
    public HashMap<Coordinate, GeomPlanarGraphEdge> edgeCoordToEdges = 
        new HashMap<Coordinate, GeomPlanarGraphEdge>();
    public Vector<Coordinate> idsToBuildings = new Vector<Coordinate>();
    public HashMap<Coordinate, String> buildingsToTypes = 
            new HashMap<Coordinate, String>();
    
    public GeomVectorField infected_agents = new GeomVectorField();
    public GeomVectorField healthy_agents = new GeomVectorField();
    ArrayList<Agent> agentList = new ArrayList<Agent>();
    
    boolean read_full_map = false; // if true, read the whole map of Seoul.
    double infect_probability = 0.5;  
    
    // variables related to distribution
    
//    ArrayList<Double> xAxisMobility = new ArrayList<>();
//    ArrayList<Double> yAxisMobility = new ArrayList<>();
    // general mobility distribution, not used currently
    
    ArrayList<Double> xAxisDistance = new ArrayList<>();
    ArrayList<Double> yAxisDistance = new ArrayList<>();
    
    ArrayList<Double> xAxisLowSpreader = new ArrayList<>();
    ArrayList<Double> yAxisLowSpreader = new ArrayList<>();
    
    ArrayList<Double> xAxisSuperSpreader = new ArrayList<>();
    ArrayList<Double> yAxisSuperSpreader = new ArrayList<>();
    
    ArrayList<Double> xAxisSymptomMove = new ArrayList<>();
    ArrayList<Double> yAxisSymptomMove = new ArrayList<>();
    
    
    /**
     * Constructor
     */
    public CovidSimulationSymptom(long seed)
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
        	// Read population and infect rate
        	Scanner s = new Scanner(new File("data/parameters.txt"));
        	population = s.nextInt();
        	infect_probability = s.nextDouble();

        	//	population = 300; // for testing
        	
        	System.out.println("rate_infection= " + infect_probability);
        	System.out.println("pop_infection= " + population);
        	s.close();
        	
            // Read the X and Y axis of the distribution from external files
        	System.out.println("reading distributions...");
        	
        	// general mobility distribution, not used currently
//            Scanner sMobility = new Scanner(new File("data/distr/xAxisMobility.txt"));
//            xAxisMobility = new ArrayList<>();
//            while(sMobility.hasNext()) {
//            	xAxisMobility.add(sMobility.nextDouble());
//            }
//            sMobility.close();
            
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
            
         // general mobility distribution, not used currently
//            sMobility = new Scanner(new File("data/distr/yAxisMobility.txt"));
//            yAxisMobility = new ArrayList<>();
//            while(sMobility.hasNext()) {
//            	yAxisMobility.add(sMobility.nextDouble());
//            }
//            sMobility.close();
            
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
        	
            // read in the roads to create the transit network
            System.out.println("reading roads layer...");

//            URL roadsFile =  Paths.get("data","road_infection_demo.shp").toUri().toURL() ; 
//            URL roadsDB =  Paths.get("data","road_infection_demo.dbf").toUri().toURL() ; 
//            URL roadsFile =  Paths.get("data","road_demo.shp").toUri().toURL() ; 
//            URL roadsDB =  Paths.get("data","road_demo.dbf").toUri().toURL() ;
            URL roadsFile =  Paths.get("data","Gangnamgu_road.shp").toUri().toURL() ; 
            URL roadsDB =  Paths.get("data","Gangnamgu_road.dbf").toUri().toURL() ;
            if (read_full_map == true)
            {
            	roadsFile =  Paths.get("data","seoul_roads_split.shp").toUri().toURL() ; 
            	roadsDB =  Paths.get("data","seoul_roads_split.dbf").toUri().toURL() ; 
            }
            ShapeFileImporter.read(roadsFile, roadsDB, roads);

            Envelope MBR = roads.getMBR();

            // read in the tracts to create the background
            System.out.println("reading tracts layer...");

//            URL areasFile = Paths.get("data","building_infection_demo.shp").toUri().toURL() ; 
//            URL areasDB = Paths.get("data","building_infection_demo.dbf").toUri().toURL() ;
//            URL areasFile = Paths.get("data","building_demo.shp").toUri().toURL() ; 
//            URL areasDB = Paths.get("data","building_demo.dbf").toUri().toURL() ;
            
            URL areasFile = Paths.get("data","Gangnamgu_building.shp").toUri().toURL() ; 
            URL areasDB = Paths.get("data","Gangnamgu_building.dbf").toUri().toURL() ;
            
            if (read_full_map == true)
            {
                areasFile = Paths.get("data","seoul_building.shp").toUri().toURL() ; 
                areasDB = Paths.get("data","seoul_building.dbf").toUri().toURL() ;            	
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
                	if (cycle_count > 14640) // 60 days.
                	{
                    	System.out.println("<log_fin>cycle="+ cycle_count+", simulated 60 days, exit.");
                    	System.exit(0);
                	}
                	
                	CovidSimulationSymptom gstate = (CovidSimulationSymptom) state;

                	// After each cycle, get all agents' location, and compare their destination building.
                	HashMap<Coordinate, ArrayList<Agent>> buildingToAgents =
                	        new HashMap<Coordinate, ArrayList<Agent>>();

                	// Agents to be isolated.
                	ArrayList<Agent> toBeRemoved = new ArrayList<Agent>();
                	
                    for (Agent a : gstate.agentList)
                    {
                    	if (gstate.agentList.size() == 0)
                    	// everyone is isolated.
                    	{
                    		System.out.println("<log_fin>Agent list is empty. All patients are isolated. Simulation finishes!");
                        	System.exit(0);
                    	}
                    	
                    	if (a.symptom == true)
                    		// if symptom is true, we don't have to check symptom_countdown. So that's why the next if needs else.
                    	{
                    		a.symptom_move_countdown --;
                    		if (a.symptom_move_countdown <= 0)
                    		{
                    			total_isolated++;
                    			System.out.println("Agent=" + a + " isolated, num_isolation=" + total_isolated);
                    			
                    			toBeRemoved.add(a);
                    			continue;
                    		}
                    		
                    	}
                    	else if (a.symptom_countdown > 0)
                    	{
                    		a.symptom_countdown --;
                    		if (a.symptom_countdown <= 0)
                    		{
                    			a.symptom = true;
                    			
                    			a.set_symptom_move_countdown(xAxisSymptomMove,yAxisSymptomMove); 
//                    			a.setMobilityAfterInfection();
                    		}
                    	}
                    	
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
                    
                    // remove isolated agents.
                    for (Agent a : toBeRemoved)
                    {
                    	gstate.agentList.remove(a);
                    	infected_agents.removeGeometry(a.getGeometry());
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
                    
                    // if everyone is infected, then stop simulation.
                    if (infectCount == population)
                    {
                    	System.out.println("<log_fin>Everyone is infected. Simulation finishes!");
                    	System.exit(0);
                    }
                    
                    // if there is no more infected case in active population, we stop simulation.
                    if (infectCount == 0)
                    {
                    	System.out.println("<log_fin>All patients are isolated. Simulation finishes!");
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
            Logger.getLogger(CovidSimulationSymptom.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex)
        {
            Logger.getLogger(CovidSimulationSymptom.class.getName()).log(Level.SEVERE, null, ex);
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

        network.createFromGeomField(roads);
        int countedge = 0;
        for (Object o : network.getEdges())
        {
            GeomPlanarGraphEdge e = (GeomPlanarGraphEdge) o;
            idsToEdges.put(countedge, e);
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

//        System.out.println("edge count " + total_edges);
        
        int count_building = 0;
        // create a hashmap, from building's Coordinate to edges(roads).
        for (Object o : buildings.getGeometries())
        {        	
            MasonGeometry b = (MasonGeometry) o;
            String type = b.getStringAttribute("type");
            Coordinate cur_b = b.getGeometry().getCoordinate();
            
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
            	System.out.println("?????");
            	System.exit(0);
            }
            
            buildingsToEdges.put(cur_b, cur_edge);
            buildingsToTypes.put(cur_b, type);
            idsToBuildings.add(cur_b);
        }
        total_buildings = buildingsToEdges.size();
//        System.exit(0);
        System.out.println("create network DONE in total" +buildingsToEdges.size() +" buildings");
    }
    
    /**
     * create an appropriate pop
     * <p/>     
     */
    public void populate()
    {
        try
        {
            for (int i = 0; i < population; i++) 
            {
            	System.out.println("in population "+i);
            	// randomly select an edge as the start position.
                int startPos = ThreadLocalRandom.current().nextInt(total_edges);
                GeomPlanarGraphEdge startingEdge =
                    idsToEdges.get(startPos);  
                 
                Agent a = new Agent(this, startingEdge);
                a.init = true;
                // super spreaders and low spreaders + mobility
                if (ThreadLocalRandom.current().nextDouble() < 0.035897 ) // 3.5897% is superspreader
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
                	a.move_speed = 0.0035972;
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
                
                // set mobility.
                // using general mobility. - currently not used.
//                a.setMobility(xAxisMobility, yAxisMobility);

                // used by Agent.setMobilityAfterInfection(). Currently not used.
//                a.xAxisLowSpreader = xAxisLowSpreader;
//                a.yAxisLowSpreader = yAxisLowSpreader;
//                
//                a.xAxisSuperSpreader = xAxisSuperSpreader;
//                a.yAxisSuperSpreader = yAxisSuperSpreader;
                
                while (!a.start(this))
                	// fail to start
                	// it is possible that agent is put on some road that is not reachable, so we try another init position.
                	// things we do in this loop is similar to new Agent().
                {
                    startPos = ThreadLocalRandom.current().nextInt(total_edges);
                    startingEdge =
                        idsToEdges.get(startPos);  
                    a.startNode = startingEdge.getDirEdge(0).getFromNode(); 

                    GeomPlanarGraphEdge goalEdge = a.chooseDestination();
                    a.endNode = goalEdge.getDirEdge(0).getToNode();
                    
                    Coordinate startCoord = null;
                    startCoord = a.startNode.getCoordinate();
                    a.updatePosition(startCoord);
                }

                MasonGeometry newGeometry = a.getGeometry();
                newGeometry.isMovable = true;
                
                // set infection
//                if (i < init_infection_rate * population)
                if (i < 1)
//                if (i < 90)
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
             }
        } catch (Exception e)
        {
            System.out.println("error in populate()");
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
    public double getCoordDistance(Coordinate a, Coordinate b)
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

