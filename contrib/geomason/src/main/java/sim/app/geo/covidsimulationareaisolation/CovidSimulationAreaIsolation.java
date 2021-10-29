/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "LICENSE" for more information
 *
 * $Id$
*/
package sim.app.geo.covidsimulationareaisolation;

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
import java.util.Random;
import java.util.Map.Entry;
import java.util.Set;
import java.util.Vector;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;

import sim.app.geo.covidsimulationareaisolation.Agent;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.field.geo.GeomVectorField;
import sim.io.geo.ShapeFileImporter;
import sim.util.geo.GeomPlanarGraph;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import java.nio.file.Files;
import java.nio.file.Path;
import org.apache.commons.math3.distribution.LogNormalDistribution;

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
public class CovidSimulationAreaIsolation extends SimState {
    private static final long serialVersionUID = 1L;
    int population = 50; // total population initialization

    /**
     * Main function allows simulation to be run in stand-alone, non-GUI mode
     */
    public static void main(String[] args) {
        doLoop(CovidSimulationAreaIsolation.class, args);
        System.exit(0);
    }
    
    public boolean init = true;
    
    public int cycle_count = 0; // used for printing log
    public int total_buildings = 0;
    public int init_infect = 1;
    public GeomVectorField buildings = new GeomVectorField();
    public int max_cycle = 14440*2;
    int total_isolated = 0; // agents that are isolated

    public boolean mitigation_switch = true;
    public int mitigation_threshold = 125;
    public boolean orig_params = true;
  
    public Vector<Coordinate> idsToBuildings = new Vector<Coordinate>();
    public HashMap<Coordinate, Integer> buildingsToVisits = 
            new HashMap<Coordinate, Integer>();
    
    // mappings for areas
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
    public HashMap<String, ArrayList<Coordinate>> typesToBuildings = 
    		new HashMap<String, ArrayList<Coordinate>>();
    // type probabilities
    public HashMap<String, HashMap<String,Double>> probType = 
    		new HashMap<String, HashMap<String,Double>>();
    public ArrayList<String> typeNames = new ArrayList<>();
    public HashMap<Integer,HashMap<String,ArrayList<Coordinate>>> areaToTypeToBldg = 
    		new HashMap<Integer,HashMap<String,ArrayList<Coordinate>>>();

    
    // agents
    public GeomVectorField infected_agents = new GeomVectorField();
    public GeomVectorField symptomatic_agents = new GeomVectorField();
    public GeomVectorField healthy_agents = new GeomVectorField();
    ArrayList<Agent> agentList = new ArrayList<Agent>();
    ArrayList<ArrayList<Agent>> symptom_cycle = new ArrayList<ArrayList<Agent>>();
    
    //params for lognormal distribution
    Double mu = 1.63; 
    Double sigma = 0.5;
    
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
    public CovidSimulationAreaIsolation(long seed) {
        super(seed);
    }

    /**
     * Initialization
     */
    @Override
    public void start() {
        super.start();
        try {
        	System.out.println("Beginning log.");
        	// Read population and infect rate
        	Scanner s = new Scanner(new File("data/parameters.txt"));
        	population = s.nextInt();
        	infect_probability = s.nextDouble();
        	caution_factor = s.nextDouble();
        	init_infect = s.nextInt();
        	mitigation_threshold=s.nextInt();
        	s.close();
        	
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
            
            // read type probability files
            
            Scanner typeNamesFile = new Scanner(new File("data/type_names.txt"));
            while (typeNamesFile.hasNext()) {
            	typeNames.add(typeNamesFile.next());
            }
            typeNamesFile.close();
            Scanner typeProbFile = new Scanner(new File("data/type_probs.txt"));
            Double prob;
            String cur_type;
            String dest_type;
        	for (int i=0; i<typeNames.size();i++) {
        		cur_type= typeNames.get(i);
        		probType.put(cur_type, new HashMap<String,Double>());
        		for (int j=0; j<typeNames.size();j++) {
        			dest_type=typeNames.get(j);
        			prob=typeProbFile.nextDouble();
                	probType.get(cur_type).put(dest_type,prob);
        		}
        		typesToBuildings.put(cur_type,new ArrayList<Coordinate>());
        	}
            typeProbFile.close();
            
            for  (int cycle_num=0;cycle_num<=max_cycle+288*300;cycle_num++) {
            	symptom_cycle.add(new ArrayList<Agent>());
            }
            
            
            // read in the tracts to create the background
            System.out.println("reading tracts layer...");
            URL areasFile = Paths.get("data","gangnam_seocho_new","gangnam_seocho_buildings.shp").toUri().toURL() ; 
            URL areasDB = Paths.get("data","gangnam_seocho_new","gangnam_seocho_buildings.dbf").toUri().toURL() ;           
            if (read_full_map == true) {
                areasFile = Paths.get("data","seoul_building_area.shp").toUri().toURL() ; 
                areasDB = Paths.get("data","seoul_building_area.dbf").toUri().toURL() ;            	
            }
            ShapeFileImporter.read(areasFile, areasDB, buildings);
            Envelope MBR = buildings.getMBR();
            createNetwork();

            // update so that everyone knows what the standard MBR is
            buildings.setMBR(MBR);

            // initialize agents
            populate();
            infected_agents.setMBR(MBR);
            healthy_agents.setMBR(MBR);
            symptomatic_agents.setMBR(MBR);
            
            // Ensure that the spatial index is updated after all the agents
            // move
            schedule.scheduleRepeating(infected_agents.scheduleSpatialIndexUpdater(), Integer.MAX_VALUE, 1.0);
            schedule.scheduleRepeating(healthy_agents.scheduleSpatialIndexUpdater(), Integer.MAX_VALUE, 1.0);
            schedule.scheduleRepeating(symptomatic_agents.scheduleSpatialIndexUpdater(), Integer.MAX_VALUE, 1.0);

            /**
             * Steppable that flips Agent paths once everyone reaches their
             * destinations
             */
            Steppable flipper = new Steppable()
            {
                @Override
                public void step(SimState state) {
                	// status to sync after each cycle.
                	
                	//make agents symptomatic
					for (Agent a: symptom_cycle.get(cycle_count)) {
						a.setSymptomatic();
						infected_agents.removeGeometry(a.getGeometry());
						symptomatic_agents.addGeometry(a.getGeometry());
						//set how long before isolation
						a.set_symptom_move_countdown(xAxisSymptomMove,yAxisSymptomMove); 
						//a.setMobilityAfterInfection();
					}
                	
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
					
                	
                	if (cycle_count > max_cycle) {
                    	System.out.println("<log_fin>cycle="+ cycle_count+", simulated"+(cycle_count/288)+"days, exit.");
                    	System.exit(0);
                	}
                	
                	CovidSimulationAreaIsolation gstate = (CovidSimulationAreaIsolation) state;

                	// After each cycle, get all agents' location, and compare their destination building.
                	HashMap<Coordinate, ArrayList<Agent>> buildingToAgents =
                	        new HashMap<Coordinate, ArrayList<Agent>>();
                	// Agents to be isolated.
                	ArrayList<Agent> toBeRemoved = new ArrayList<Agent>();
                	
                    for (Agent a : gstate.agentList)
                    {         
                    	if (gstate.agentList.size() == 0) { // everyone is isolated.
                        		System.out.println("<log_fin>Agent list is empty. All patients are isolated. Simulation finishes!");
                            	System.exit(0);
                        	}
                    	if (a.checkSymptomatic()) {
                    		a.symptom_move_countdown --;
                    		if (a.symptom_move_countdown <= 0) {
                    			total_isolated++;
                    			System.out.println("Agent=" + a + " isolated, num_isolation=" + total_isolated);
                    			toBeRemoved.add(a);
                    			continue;
                    		}	
                    	}
                    	
                    	if ((a.getInBuildingStatus() == true) && (a.atHome == false)) { // a is inside a building
                    		Coordinate curBuilding = a.getDestinationBuilding();
                    		ArrayList<Agent> tmp = buildingToAgents.get(curBuilding);
                    		// check the status of this building: is there any agent also in this building?
                    		if (tmp == null) {
                    			tmp = new ArrayList<Agent>();
                    		}
                    		tmp.add(a);
                    		buildingToAgents.put(curBuilding, tmp);
                    		// print log: we find a in the building.
                    		System.out.println("<log> cycle="+ cycle_count + " agent="+ a + " building="+ curBuilding + " infection=" + a.checkInfection());
                    	} else if (a.atHome == true) {
                    		System.out.println("<log> cycle="+ cycle_count + " agent="+ a + " home building="+ a.homeBldg + " infection=" + a.checkInfection());
                    	}
                    	
                    }
                    
                    // remove isolated agents.
                    for (Agent a : toBeRemoved) {
                    	gstate.agentList.remove(a);
                    	symptomatic_agents.removeGeometry(a.getGeometry());
                    	System.out.println("Agent "+a+" isolated");
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
                    for (Agent a : gstate.agentList) {
                    	if (a.checkInfection() == true) {
                    		infectCount ++;
                    	}
                    }
                    System.out.println("<log> cycle="+ cycle_count + " infection="+ infectCount + " isolation="+total_isolated + " total cases="+(infectCount+total_isolated));
                    
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
                    if (infectCount == population) {
                    	System.out.println("<log_fin>Everyone is infected. Simulation finishes!");
                    	System.exit(0);
                    }
                    // if there is no more infected case in active population, we stop simulation.
                    if (infectCount == 0) {
                    	System.out.println("<log_fin>All patients are isolated. Simulation finishes!");
                    	System.exit(0);
                    }
                }
            };
            schedule.scheduleRepeating(flipper, 1);

        } catch (FileNotFoundException ex) {
            System.out.println("Error: missing required data file");
        } catch (IOException ex) {
            Logger.getLogger(CovidSimulationAreaIsolation.class.getName()).log(Level.SEVERE, null, ex);
        } catch (Exception ex) {
            Logger.getLogger(CovidSimulationAreaIsolation.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void infect(ArrayList<Agent> curAgents) {
    	int infected = 0;
    	for (Agent a : curAgents) {
   			if (a.checkSymptomatic() == true) {
   				infected ++;
    		}
    	}
    
   		double pr_infection = 1; 
   		for (int i = 0; i < infected; i++) {
    		pr_infection = (1-infect_probability) * pr_infection;
    	}
    
    	pr_infection = 1 - pr_infection;
    
    	if (infected > 0) {
    		for (Agent a : curAgents) {	
    			if (a.checkInfection() == false) {
    				double tmp = ThreadLocalRandom.current().nextDouble();
    				if (tmp < pr_infection) {
    					a.setInfection();
    					healthy_agents.removeGeometry(a.getGeometry());
    					infected_agents.addGeometry(a.getGeometry());
    					
    					LogNormalDistribution logNormal = new LogNormalDistribution(sigma, mu);
    				    Double num_days = logNormal.sample();
    				    while ((num_days>20.0) || (num_days<1.0)) {
    				    	num_days = logNormal.sample();
    				    }
    				    int num_cycles_sympt = (int) Math.ceil(num_days * 288);
    				    
    					symptom_cycle.get(cycle_count+num_cycles_sympt).add(a);
    					System.out.println("<log-infection> cycle="+ cycle_count + " agent="+ a);
    				}
    			}
   			}
   		}
    }

    /**
     * Create the network the agents will traverse
     * <p/>
     */
    private void createNetwork() {
    	init = true;
        System.out.println("creating network...");
        
        int a_count=0;
    	ArrayList<Double> origin_areas = new ArrayList<Double>();
        Scanner origin_probs_file;
		try {
			origin_probs_file = new Scanner(new File("data/distr/origin_probs.txt")); //just to get number of areas...will replace
			while(origin_probs_file.hasNext()) {
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
    		//areasToEdges.put(i,new ArrayList<GeomPlanarGraphEdge>());
    		areasToHotspots.put(i, new ArrayList<Coordinate>());
    		areaToTypeToBldg.put(i, new HashMap<String,ArrayList<Coordinate>>());
    	}
    	HashMap<String,String> shptypeToType = new HashMap<String,String>();
        HashMap<Integer,ArrayList<Double>> miscTypeProbs = new HashMap<Integer,ArrayList<Double>>();

    	try {
    		Scanner shptypeFile = new Scanner(new File("data/shptype_to_type.txt"));
    		while (shptypeFile.hasNext()) {
    			String values = shptypeFile.next();
    			shptypeToType.put(values.split(",")[0], values.split(",")[1]);
    		}
    		shptypeFile.close();
    		Double prob=0.0;
            Double[] prob_sum= {0.0,0.0};
        	Scanner typeMiscFile = new Scanner(new File("data/misc_types.txt"));
            for (int i=0;i<typeNames.size();i++) {
            	String from_file = typeMiscFile.next();
            	String[] split_data = from_file.split(",");
            	//String cur_misc_type = typeNames.get(i);
            	for (int j=0;j<2;j++) {
            		prob = Double.valueOf(split_data[j]);
            		prob_sum[j]+=prob;
            		if (!miscTypeProbs.containsKey(j)) {
            			miscTypeProbs.put(j,new ArrayList<Double>());
            		}
            		miscTypeProbs.get(j).add(prob_sum[j]);
            	}
            }
            typeMiscFile.close();
    	} catch (FileNotFoundException e1) {
    		System.out.println("No mapping file for shapefile types to types");
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
		int misc_count=0;
        for (Object o : buildings.getGeometries()) {        	
            MasonGeometry b = (MasonGeometry) o;
            String shptype = b.getStringAttribute("type"); 
            Coordinate cur_b = b.getGeometry().getCoordinate();
          
            int area_code=b.getIntegerAttribute("area_code");
            String type;
            if (shptypeToType.containsKey(shptype)) {
            	type = shptypeToType.get(shptype); 
            } else {
            	Random rndUnif = new Random();
            	double u = rndUnif.nextDouble();
            	int idx = 0;
            	Double cur_prob = 0.0;
            	String misc_type = typeNames.get(idx);
            	while (u > cur_prob) {
            		try {
	            		cur_prob = miscTypeProbs.get(area_code).get(idx);
	            		misc_type = typeNames.get(idx);
	            		idx+=1; 
            		} catch (Exception e) {
            			misc_type = typeNames.get(idx-1);
            			break;
            		}
            	}
            	type = misc_type;
            	misc_count++;
            }
            if (!areaToTypeToBldg.get(area_code).containsKey(type)) {
            	areaToTypeToBldg.get(area_code).put(type,new ArrayList<Coordinate>());
            }
            areaToTypeToBldg.get(area_code).get(type).add(cur_b);
            buildingsToVisits.put(cur_b, 0);
            buildingsToAreas.put(cur_b,area_code);
            buildingsToTypes.put(cur_b, type);
            typesToBuildings.get(type).add(cur_b);
            idsToBuildings.add(cur_b);
            areasToBuildings.get(area_code).add(cur_b);
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

        total_buildings = idsToBuildings.size();
        System.out.println("create network DONE in total " +total_buildings +" buildings");
        //System.exit(0);
        /*System.out.println("-----------------");

        System.out.println(areaToTypeToBldg.get(0).get("work_rec").size());
        System.out.println(areaToTypeToBldg.get(0).get("rec_night").size());
        System.out.println(areaToTypeToBldg.get(0).get("residence").size());
        System.out.println(areaToTypeToBldg.get(0).get("medical").size());
        System.out.println(areaToTypeToBldg.get(0).get("worship").size());
        System.out.println(areaToTypeToBldg.get(0).get("work_essential").size());
        System.out.println("-----------------");

        System.out.println(areaToTypeToBldg.get(1).get("work_rec").size());
        System.out.println(areaToTypeToBldg.get(1).get("rec_night").size());
        System.out.println(areaToTypeToBldg.get(1).get("residence").size());
        System.out.println(areaToTypeToBldg.get(1).get("medical").size());
        System.out.println(areaToTypeToBldg.get(1).get("worship").size());
        System.out.println(areaToTypeToBldg.get(1).get("work_essential").size());
        System.out.println("-----------------");
        System.out.println(probType);
        System.out.println("-----------------");
        System.out.println(misc_count);
        System.exit(0);*/
    }
    
    /**
     * create an appropriate pop
     * <p/>     
     */
    public void populate() {
        try {
        	Scanner origin_probs_file= new Scanner(new File("data/distr/origin_probs.txt"));
        	ArrayList<Double> origin_probs = new ArrayList<Double>();
        	int a_count=0;
        	ArrayList<Double> origin_areas = new ArrayList<Double>();
        	while(origin_probs_file.hasNext()) {
        		origin_probs.add(origin_probs_file.nextDouble());
        		origin_areas.add(Double.valueOf(a_count));
        		a_count++;
        	}
        	origin_probs_file.close();
        	Scanner origin_superspreader = new Scanner(new File("data/area_superspreader.txt"));
        	ArrayList<Double> superspreader_probs = new ArrayList<Double>();
        	while(origin_superspreader.hasNext()) {
        		superspreader_probs.add(origin_superspreader.nextDouble());
        	}
        	int num_areas=a_count;
        	Scanner area_probs = new Scanner(new File("data/distr/area_probs.txt"));
        	ArrayList<ArrayList<Double>> area_probs_all = new ArrayList<ArrayList<Double>>();
        	
        	for (int j=0;j<num_areas;j++) {
        		area_probs_all.add(new ArrayList<Double>());
        		for (int k=0;k<num_areas;k++) {
        			area_probs_all.get(j).add(area_probs.nextDouble());
        		}
        	}
        	area_probs.close();

        	//String fileContent="";
        	for (int i = 0; i < population; i++) 
            {
            	System.out.println("in population "+i);
            	Agent a = new Agent(this);
                a.init = true;
                a.setArea(origin_areas, origin_probs,  area_probs_all);
                // super spreaders and low spreaders + mobility
                if (ThreadLocalRandom.current().nextDouble() < superspreader_probs.get(a.area) ) {
                	a.superspreader = true;
                	a.setMobility(xAxisSuperSpreader, yAxisSuperSpreader);
                	System.out.println("Agent="+a+" set superspreader");
                } else {
                	a.setMobility(xAxisLowSpreader, yAxisLowSpreader);
                	System.out.println("Agent="+a+" set lowspreader");
                }
                
                // speed
                if (i % 2 == 0) {
                	// walk
                	a.speed = 3.0/12.0;
                	a.move_speed = 3.0/12.0;
                } else {
                	// car
                	//10-25
                	double maxspeed = 25;
                	double minspeed = 10;
                	double tmp = ThreadLocalRandom.current().nextDouble();
                	a.speed = (minspeed + (maxspeed-minspeed)*tmp)/12;
                	a.move_speed = a.speed;
                }
                ArrayList<Coordinate> blist = areasToBuildings.get(a.checkArea());
               	int b_pos=ThreadLocalRandom.current().nextInt(blist.size());
               	a.homeBldg=blist.get(b_pos);
               	a.destination_building=a.homeBldg;
               	a.prev_building=a.homeBldg;
               	a.updatePosition(a.homeBldg);
               	//a.destination_building=a.chooseDestination();
               	while(!a.start(this)) {
               		b_pos=ThreadLocalRandom.current().nextInt(blist.size());
               		a.homeBldg=blist.get(b_pos);
               		a.destination_building=a.homeBldg;
               		a.prev_building=a.homeBldg;
               		a.updatePosition(a.homeBldg);
               		//a.destination_building=a.chooseDestination();
               	}
                MasonGeometry newGeometry = a.getGeometry();
                newGeometry.isMovable = true;
                
                // set infection
                if (i < init_infect) {
                	a.setInfection();
                	
                	/*LogNormalDistribution logNormal = new LogNormalDistribution(sigma, mu);
				    Double num_days = logNormal.sample();
				    while (num_days>20.0) {
				    	num_days = logNormal.sample();
				    }
				    int num_cycles_sympt = (int) Math.ceil(num_days * 288);*/
                	int num_cycles_sympt = 1; //make contagious next cycle
				    
					symptom_cycle.get(num_cycles_sympt).add(a);
                	infected_agents.addGeometry(newGeometry);
                } else {
                	healthy_agents.addGeometry(newGeometry);
                }
 
                agentList.add(a);
                schedule.scheduleRepeating(a);
                
                init = false;
                a.init = false;
                //fileContent=fileContent+String.valueOf(a.getGeometry())+"\t";
                //fileContent=fileContent+String.valueOf(a.superspreader)+"\t";
                //fileContent=fileContent+String.valueOf(a.getMobility())+"\t";
                //fileContent=fileContent+String.valueOf(a.checkInfection())+"\t";
                //fileContent=fileContent+String.valueOf(a.checkArea())+"\t";
                //fileContent=fileContent+String.valueOf(a.move_speed)+"\n";
                
                //Path fileName = Path.of("data/area_profiling_pop10000.txt");                
                //Files.writeString(fileName, fileContent);
                
             }
        } catch (Exception e) {
        	System.out.println(e);
            System.out.println("error in populate()");
            System.exit(1);
        }
    }
 
    public double getCoordDistance(Coordinate a, Coordinate b) { //returns distance in meters
    	double lng1 = a.x;    // longitude
    	double lng2 = b.x;
    	double lat1 = a.y;
    	double lat2 = b.y;
    
    	double x = (lng2 - lng1) * Math.PI * 6.371229 * 1e6 * Math.cos(((lat1 + lat2) / 2) * Math.PI / 180) / 180;
        double y = (lat1 - lat2) * Math.PI * 6.371229 * 1e6 / 180;
        return Math.hypot(x, y);
    }
}

