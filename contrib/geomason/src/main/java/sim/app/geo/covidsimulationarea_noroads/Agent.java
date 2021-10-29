/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "LICENSE" for more information
 *
 * $Id$
 */
package sim.app.geo.covidsimulationarea_noroads;
import java.lang.Math;
import com.vividsolutions.jts.geom.*;
import com.vividsolutions.jts.linearref.LengthIndexedLine;
import com.vividsolutions.jts.planargraph.Node;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.Map.Entry;
import java.util.concurrent.ThreadLocalRandom;

import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.geo.GeomPlanarGraphDirectedEdge;
import sim.util.geo.GeomPlanarGraphEdge;
import sim.util.geo.MasonGeometry;
import sim.util.geo.PointMoveTo;

/**
 *  Our simple agent for the CampusWorld GeoMASON example.  The agent randomly wanders
 *   around the campus walkways.  When
 *  the agent reaches an intersection, it chooses a random direction and continues on.   
 *
 */
@SuppressWarnings("restriction")
public final class Agent implements Steppable
{
	public boolean init = true;
	public boolean superspreader = false; // true: superspreader; false: low spreader
    private boolean infected = false; // flag of whether this agent has covid.
    public int area = -1;
    public ArrayList<Double> area_probs= new ArrayList<Double>();
    public ArrayList<Double> area_list= new ArrayList<Double>();
    private double move_probability = 0.9; // probability of leaving a building when an agent is inside the building.
    private static final long serialVersionUID = -1113018274619047013L;
    CovidSimulationAreaNoRoads world;    
    // Residence/Work Attributes
    public Coordinate destination_building = null;
    public Coordinate prev_building = null;   
    public Coordinate homeBldg = null;
    // point that denotes agent's position   
    private MasonGeometry location;
    boolean inBuilding = false;
    public boolean atHome = true;
    PointMoveTo pointMoveTo = new PointMoveTo();
    public double move_speed = 0.0005;  // move_speed is a fixed number of speed
    public double speed = move_speed; // speed is the actual moved distance in each cycle (direction + tail if too close to a building)
    public int stay_cycle = 3; // Agents have to stay in a building for at least 15 minutes.
    int cur_stay_cycle = stay_cycle;
    int mv_cycle = 0;
    /** This is the wrapper object in the agents layer.  We need a handle on
     * it so that we can update our location with each step().
     */

    /** Constructor Function */
    public Agent(CovidSimulationAreaNoRoads g) {
        world = g;
        GeometryFactory fact = new GeometryFactory();
        location = new MasonGeometry(fact.createPoint(new Coordinate(10, 10))) ;
    }

    /** Initialization of an Agent: find an A* path from startNode to endNode
     *
     * @param state
     * @return whether or not the agent successfully found a path
     */
    public boolean start(CovidSimulationAreaNoRoads state) {
    	//aks what cases might we need to select new destinations?    
    	return true;
    }

    /** Called every tick by the scheduler */
    /** moves the agent along the path */
    public void step(SimState state) {
    	if (inBuilding == true) {
    		if (cur_stay_cycle > 0) {
    			cur_stay_cycle -= 1;
    		} else {
    			double cur_rand = ThreadLocalRandom.current().nextDouble();
                double mv_prob;
                if (destination_building == homeBldg) {
                	mv_prob=move_probability/world.caution_factor;
                }else {
                	mv_prob=move_probability;
                }
                if (cur_rand < mv_prob) { // leaving the building.
                	inBuilding = false;
                	atHome = false;
                	prev_building = destination_building;
                	destination_building = chooseDestination();
                    
                    while (!this.start(world)) { // start is necessary, to update the path.
                        world.buildingsToVisits.put(destination_building, Integer.valueOf(world.buildingsToVisits.get(destination_building)-1));
                    	destination_building = chooseDestination();
                    }
                    Double dist_scale = ThreadLocalRandom.current().nextDouble()+1.0;
                    mv_cycle = (int) (0.000621371 * Math.ceil((getCoordDistance(prev_building,destination_building)*dist_scale) / speed));
                }    			
    		}
    	} else {
    		if (mv_cycle > 0) {
    			mv_cycle -= 1;
    		} else {
    			//we have arrived
    			inBuilding = true;
    			if (destination_building == homeBldg) {
    				atHome = true;
    			} else {
    				atHome = false;
    			}
    			stay_cycle=3;
    			cur_stay_cycle = stay_cycle;
    			updatePosition(destination_building);
    		}
    	}
    }                
 
    ///////////// HELPER FUNCTIONS //////////////////////////// 
    public boolean checkInfection() {
        return infected;
    }    
    public void setInfection() { // for return value: true - success; false - failed
		System.out.println("in Agent: " +  this + " infected.");
    	infected = true;
    }
   
    public void setArea(ArrayList<Double> xAxis, ArrayList<Double> yAxis, ArrayList<ArrayList<Double>> mv_distr) {
    	Double res = sampleFromDistr(xAxis,yAxis);
    	area = res.intValue();
    	area_list=xAxis;
    	area_probs = mv_distr.get(area);  	
    }
    public int checkArea() {
        return this.area;
    }  
    public ArrayList<Double> checkAreaProbs() {
    	return this.area_probs;
    }
  
    public void setMobility (ArrayList<Double> xAxis, ArrayList<Double> yAxis) {
    	double mobi = sampleFromDistr(xAxis, yAxis)+0.001;
    	move_probability = mobi;
    }
    public double getMobility() {
    	return this.move_probability;
    }
    public double getSpeed() {
    	return this.speed;
    }
    /** move the agent to the given coordinates */
    public void updatePosition(Coordinate c) {
        pointMoveTo.setCoordinate(c);
        if (infected == true) {
        	world.infected_agents.setGeometryLocation(location, pointMoveTo);
        } else {
        	world.healthy_agents.setGeometryLocation(location, pointMoveTo);
        }
        
        if (location.geometry.getCoordinate() != c) {
            location.geometry.apply(pointMoveTo);     	
        }   
    }

    /** return geometry representing agent location */
    public MasonGeometry getGeometry() {
        return location;
    }
    
    public double getCoordDistance(Coordinate a, Coordinate b) {
    	double lng1 = a.x;    // longitude
    	double lng2 = b.x;
    	double lat1 = a.y;
    	double lat2 = b.y;
    
    	double x = (lng2 - lng1) * Math.PI * 6.371229 * 1e6 * Math.cos(((lat1 + lat2) / 2) * Math.PI / 180) / 180;
        double y = (lat1 - lat2) * Math.PI * 6.371229 * 1e6 / 180;
        return Math.hypot(x, y);
    }
    
    public Coordinate chooseDestination() { // return the id of the destination edge.
    	Coordinate nextBuilding = new Coordinate();
    	//determine if return home
    	Random rand_num = new Random();
    	double r = rand_num.nextDouble();
    	if (r>move_probability) {
    		//return to home node
    		nextBuilding=homeBldg;
    	} 
    	else {
    		ArrayList<Double> enumTypeX = new ArrayList<>();
    		ArrayList<Double> typeProbY = new ArrayList<>();
    		for (int i=0;i<world.typeNames.size();i++) {
    			enumTypeX.add(Double.valueOf(i));
    			String cur_type = world.buildingsToTypes.get(destination_building);
    			typeProbY.add(world.probType.get(cur_type).get(world.typeNames.get(i)));
    		}
    		String nextType = world.typeNames.get( (int) sampleFromDistr(enumTypeX,typeProbY));
    		int nextArea = (int) sampleFromDistr(area_list,area_probs);
    		
    		ArrayList<Coordinate> bldgs_type_area = world.areaToTypeToBldg.get(nextArea).get(nextType);
    		ArrayList<Coordinate> area_hotspots = world.areasToHotspots.get(nextArea);
    		double hotspot_prob=0.0;
    		ArrayList<Double> problist = new ArrayList<Double>();
    		for (int i=0;i<area_hotspots.size();i++) {
    			hotspot_prob=hotspot_prob+world.hotspotsToProbs.get(area_hotspots.get(i));
    			problist.add(hotspot_prob);
    		}
    		Random rndUnif = new Random();
    		double u = rndUnif.nextDouble();
    	
    		if (u<=hotspot_prob) {
    			Coordinate spot = area_hotspots.get(area_hotspots.size()-1);
    			for (int j=0;j<problist.size();j++) {   
    				if (u<=problist.get(j)) {
    					spot = area_hotspots.get(j);
    					break;
    				}
    			}
    			int num_bldgs = world.hotspotsToBuildings.get(spot).size();
    			int nextPos = ThreadLocalRandom.current().nextInt(num_bldgs);
    			nextBuilding = world.hotspotsToBuildings.get(spot).get(nextPos);
    		}
    		else {
    			int num_bldgs = bldgs_type_area.size();
    			int nextPos = ThreadLocalRandom.current().nextInt(num_bldgs);
    			nextBuilding = bldgs_type_area.get(nextPos);

    		}
    	}
        world.buildingsToVisits.put(nextBuilding, Integer.valueOf(world.buildingsToVisits.get(nextBuilding)+1));
        return nextBuilding;
    }
    
    //Get a random sample from a distribution with PDF defined by xAxis and yAxis
    public double sampleFromDistr(ArrayList<Double> xAxis, ArrayList<Double> yAxis) {
    	//Compute the sum of all items on x-axis
    	float total = 0;
    	for(int i=0; i<yAxis.size(); i++) {
    		total += yAxis.get(i);
    	}
    	//Normalize the x-axis over its sum to compute the probability of an event (as a percentage)
    	ArrayList<Double> probList = new ArrayList<>();
    	for(int i=0; i<yAxis.size(); i++) {
    		probList.add(i, yAxis.get(i)/total);
    	}
    	//Compute the CDF by summing (cumsum in Python) all probabilities in probList
    	ArrayList<Double> cdf = new ArrayList<>();
    	for(int i=0; i<probList.size(); i++) {
    		if(i>0) {
    			double prev = cdf.get(cdf.size()-1);
    			cdf.add(i, prev+probList.get(i));
    		} else {
    			cdf.add(i, probList.get(i));
    		}
    	}
    	
    	//Get a random number from Uniform(0,1)
    	Random rndUnif = new Random();
    	double u = rndUnif.nextDouble();
    	//Get a random number from the distribution
    	for(int i=0; i<cdf.size(); i++) {
    		if(u < cdf.get(i)) {
    			return xAxis.get(i);
    		}
    	}
    	return -1; //Return -1 if, for any reason (an error!), a random number could not be obtained from the distribution
    }

    public Coordinate getDestinationBuilding() {
    	return destination_building;
    }    
    public boolean getInBuildingStatus() {
    	return inBuilding;
    }
}

