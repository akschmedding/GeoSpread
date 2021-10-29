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
    public boolean area_stay = true;
    public ArrayList<Double> area_probs= new ArrayList<Double>();
    public ArrayList<Double> area_list= new ArrayList<Double>();
    private double move_probability = 0.9;
    // probability of leaving a building when an agent is inside the building.

    private static final long serialVersionUID = -1113018274619047013L;
    CovidSimulationArea world;
    
    // Residence/Work Attributes
    Coordinate destination_building;
    
    public Node startNode = null;
    public Node endNode = null;
    
    public Coordinate homeBldg = null;
    // point that denotes agent's position
    
    private MasonGeometry location;
    
    // How much to move the agent by in each step()
    private double moveRate = .001;
    
    // Used by agent to walk along line segment
    private LengthIndexedLine segment = null;
    double startIndex = 0.0; // start position of current line
    double endIndex = 0.0; // end position of current line
    double currentIndex = 0.0; // current location along line
    GeomPlanarGraphEdge currentEdge = null;
    int linkDirection = 1;
    GeomPlanarGraphEdge startingEdge;
    
    ArrayList<GeomPlanarGraphDirectedEdge> pathFromStartToEnd =
        new ArrayList<GeomPlanarGraphDirectedEdge>();
    int indexOnPath = 0;
    int pathDirection = 1;
    boolean reachedDestination = false;
    boolean inBuilding = false;
    PointMoveTo pointMoveTo = new PointMoveTo();
    
    public double move_speed = 0.0005;  // move_speed is a fixed number of speed
    public double speed = move_speed; // speed is the actual moved distance in each cycle (direction + tail if too close to a building)
    public int stay_cycle = 3; // Agents have to stay in a building for at least 15 minutes.
    int cur_stay_cycle = stay_cycle;
    
    /** This is the wrapper object in the agents layer.  We need a handle on
     * it so that we can update our location with each step().
     */

    /** Constructor Function */
    public Agent(CovidSimulationArea g, 
                 GeomPlanarGraphEdge startEdge)
    {
        world = g;
        startingEdge=startEdge;
        startNode = startingEdge.getDirEdge(0).getFromNode();
        
        GeometryFactory fact = new GeometryFactory();
        location = new MasonGeometry(fact.createPoint(new Coordinate(10, 10))) ;
       
    }

    /** Initialization of an Agent: find an A* path from startNode to endNode
     *
     * @param state
     * @return whether or not the agent successfully found a path
     */
    public boolean start(CovidSimulationArea state)
    {
    	inBuilding = false; // agent starts move.
    	cur_stay_cycle = stay_cycle;
        if (!findNewAStarPath(state))
        {
        	return false;
        }
        
        if (pathFromStartToEnd.isEmpty())
        {
        	System.out.println("pathFromStartToEnd.isEmpty()");
            return false;
        } 
        else
        {
            return true;
        }
    }

    /** Plots a path between the Agent's home Node and its work Node */
    private boolean findNewAStarPath(CovidSimulationArea geoTest)
    {
        // get the home and work Nodes with which this Agent is associated
        Node currentJunction = geoTest.network.findNode(location.geometry.getCoordinate());
        Node destinationJunction = endNode;

        while (currentJunction == null)
        {
        	Coordinate cur_b = location.geometry.getCoordinate();
        	GeomPlanarGraphEdge cur_edge = new GeomPlanarGraphEdge(null);
            double min_dist = 1000000000;
            
            // go through building coordinates to find the closest
            Set<Entry<Coordinate, GeomPlanarGraphEdge>> coords = world.edgeCoordToEdges.entrySet();
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
            
            Node cur_node = cur_edge.getDirEdge(0).getToNode();
            Coordinate cur_coord = cur_node.getCoordinate();
            updatePosition(cur_coord);
            currentJunction = cur_node;
        }
 
        // find the appropriate A* path between them
        AStar pathfinder = new AStar();
        ArrayList<GeomPlanarGraphDirectedEdge> path =
            pathfinder.astarPath(currentJunction, destinationJunction);

        if (geoTest.init == true || init == true)
        {
        	if (path == null || path.size() <= 0) {
        		return false;
        	}
        }
    		
        double limit_dist = 20;
        // if the path works, lay it in
        while (path == null || path.size() <= 0)
        {
        	if (limit_dist > 1000) {
        		
        		return false;
        	}
        	if (geoTest.init == true || init == true) {
        		return false;
        	}

        	Coordinate cur_b = location.geometry.getCoordinate();
        	GeomPlanarGraphEdge cur_edge = new GeomPlanarGraphEdge(null);
            double min_dist = 1000000000;
            
            limit_dist = limit_dist + 20;
            // go through building coordinates to find the closest
            Set<Entry<Coordinate, GeomPlanarGraphEdge>> coords = world.edgeCoordToEdges.entrySet();
            Iterator<Entry<Coordinate, GeomPlanarGraphEdge>> iterator = coords.iterator();
            double cur_dist=-1.0;
            while(iterator.hasNext()) {
            	Map.Entry<Coordinate, GeomPlanarGraphEdge> mentry = (Map.Entry<Coordinate, GeomPlanarGraphEdge>)iterator.next();
	               
	            if (mentry.getKey() == cur_b)
	            	continue;
	               
	            if (mentry.getKey() == currentJunction.getCoordinate())
	            	continue;
	            		   
	            cur_dist = getCoordDistance(cur_b, mentry.getKey());
	            if (cur_dist < limit_dist)
	            	continue;
	            if (cur_dist < min_dist)
	            {
	            	min_dist = cur_dist;
	            	cur_edge = mentry.getValue();
	            }
            }
            Node cur_node;
                        
            cur_node = cur_edge.getDirEdge(0).getToNode();
            Coordinate cur_coord = cur_node.getCoordinate();
            if (cur_coord == cur_b)
            {
            	 cur_node = cur_edge.getDirEdge(1).getToNode();
            	 cur_coord = cur_node.getCoordinate();
            }
            updatePosition(cur_coord);
            currentJunction = cur_node;
            path = pathfinder.astarPath(currentJunction, destinationJunction);
        }

        // save it
        pathFromStartToEnd = path;

        // set up how to traverse this first link
        GeomPlanarGraphEdge edge =
            (GeomPlanarGraphEdge) path.get(0).getEdge();
        setupEdge(edge);

        // update the current position for this link
        updatePosition(segment.extractPoint(currentIndex));
        return true;
    }


    /** Called every tick by the scheduler */
    /** moves the agent along the path */
    public void step(SimState state)
    {
        // check that we've been placed on an Edge
        if (segment == null)
        {
            System.out.println(this + " segment null"); // shouldn't happen
            return;
        } // check that we haven't already reached our destination
        else if (reachedDestination)
        {
            return;
        }

        // make sure that we're heading in the right direction
        if (pathDirection < 0)
        {
            flipPath();
        }

        // move along the current segment
        speed = linkDirection*moveRate;
        
        currentIndex += speed;

        // check to see if the progress has taken the current index beyond its goal
        // given the direction of movement. If so, proceed to the next edge
        if (linkDirection == 1 && currentIndex > endIndex)
        {
            Coordinate currentPos = segment.extractPoint(endIndex);
            updatePosition(currentPos);
            transitionToNextEdge(currentIndex - endIndex);
        } else if (linkDirection == -1 && currentIndex < startIndex)
        {
            Coordinate currentPos = segment.extractPoint(startIndex);
            updatePosition(currentPos);
            transitionToNextEdge(startIndex - currentIndex);
        } else
        { // just update the position!
            Coordinate currentPos = segment.extractPoint(currentIndex);

            updatePosition(currentPos);
        }
    }

    /** Flip the agent's path around */
    void flipPath()
    {
        reachedDestination = false;
        pathDirection = -pathDirection;
        linkDirection = -linkDirection;
    }

    /**
     * Transition to the next edge in the path
     * @param residualMove the amount of distance the agent can still travel
     * this turn
     */
    void transitionToNextEdge(double residualMove)
    {
        // update the counter for where the index on the path is
        indexOnPath += pathDirection;

        // check to make sure the Agent has not reached the end
        // of the path already
        if ((pathDirection > 0 && indexOnPath >= pathFromStartToEnd.size())
            || (pathDirection < 0 && indexOnPath < 0))// depends on where you're going!
        {
        	// Arrived.
         	inBuilding = true;
        	
        	if (cur_stay_cycle > 0)
        	{
        		// stay for at least 3 cycles
        		cur_stay_cycle --;
        		return;
        	}
        	
            double cur_rand = ThreadLocalRandom.current().nextDouble();
            double mv_prob;
            if (destination_building == homeBldg) {
            	mv_prob=move_probability/world.caution_factor;
            }else {
            	mv_prob=move_probability;
            }
            if (cur_rand < mv_prob)
            	// leaving the building.
            {
            	indexOnPath = 0;
            	
            	startNode =  endNode; // cur_pos = end_pos = endNode = next_start_pos = startNode
                GeomPlanarGraphEdge goalEdge = chooseDestination();
                endNode = goalEdge.getDirEdge(0).getToNode();

                Coordinate startCoord = startNode.getCoordinate();
                updatePosition(startCoord);
                
                while (!this.start(world)) // start is necessary, to update the path.
                {
                	// if cannot reach, choose another end node.
                    world.buildingsToVisits.put(destination_building, Integer.valueOf(world.buildingsToVisits.get(destination_building)-1));
                    //System.out.println(destination_building);
                    goalEdge = chooseDestination();
                    endNode = goalEdge.getDirEdge(0).getToNode();
                }
            }
            return; // agents start moving at the next cycle, so this cycle agents do not move.
        }

        // move to the next edge in the path
        GeomPlanarGraphEdge edge =
            (GeomPlanarGraphEdge) pathFromStartToEnd.get(indexOnPath).getEdge();
        setupEdge(edge);
        speed = linkDirection*residualMove;
        currentIndex += speed;

        // check to see if the progress has taken the current index beyond its goal
        // given the direction of movement. If so, proceed to the next edge
        if (linkDirection == 1 && currentIndex > endIndex)
        {
            transitionToNextEdge(currentIndex - endIndex);
        } else if (linkDirection == -1 && currentIndex < startIndex)
        {        	
            transitionToNextEdge(startIndex - currentIndex);
        }
    }

    ///////////// HELPER FUNCTIONS ////////////////////////////

    /** Sets the Agent up to proceed along an Edge
     * @param edge the GeomPlanarGraphEdge to traverse next
     * */
    void setupEdge(GeomPlanarGraphEdge edge)
    {
        // set up the new segment and index info
        LineString line = edge.getLine();
        segment = new LengthIndexedLine(line);
        startIndex = segment.getStartIndex();
        endIndex = segment.getEndIndex();
        linkDirection = 1;

        // check to ensure that Agent is moving in the right direction
        double distanceToStart = line.getStartPoint().distance(location.geometry),
            distanceToEnd = line.getEndPoint().distance(location.geometry);
        if (distanceToStart <= distanceToEnd)
        { // closer to start
            currentIndex = startIndex;
            linkDirection = 1;
        } else if (distanceToEnd < distanceToStart)
        { // closer to end
            currentIndex = endIndex;
            linkDirection = -1;
        }
    }

    public boolean checkInfection()
    {
        return infected;
    }
    
    public void setInfection()
    // for return value: true - success; false - failed
    {
		System.out.println("in Agent: " +  this + " infected.");

    	infected = true;
    }
    public void setArea(ArrayList<Double> xAxis, ArrayList<Double> yAxis, ArrayList<Double> stay_probs, ArrayList<ArrayList<Double>> mv_distr)
    {
    	Double res = sampleFromDistr(xAxis,yAxis);
    	area = res.intValue();
    	area_list=xAxis;
    	Random rndUnif = new Random();
    	double u = rndUnif.nextDouble();
    	//Get a random number from the distribution
    	if (u<stay_probs.get(area)) 
    	{
    		area_stay=true;
    		for (int i=0;i<xAxis.size();i++) 
    		{
    			if (i!=area) {
    				area_probs.add(0.0);
    			}
    			else {
    				area_probs.add(1.0);
    			}
    		}
    	}
    	else
    	{
    		area_stay=false;
    		area_probs = mv_distr.get(area);
    	}
    	 
    	
    }
    public int checkArea()
    {
        return this.area;
    }
   
    public ArrayList<Double> checkAreaProbs()
    {
    	return this.area_probs;
    }
    public boolean checkAreaStay() {
    	return this.area_stay;
    }
  
    public void setMobility (ArrayList<Double> xAxis, ArrayList<Double> yAxis)
    {
    	double mobi = sampleFromDistr(xAxis, yAxis)+0.001;
    	move_probability = mobi;
    }
    
    public double getMobility() {
    	return this.move_probability;
    }
    
    /** move the agent to the given coordinates */
    public void updatePosition(Coordinate c)
    {
        pointMoveTo.setCoordinate(c);
        if (infected == true)
        {
        	world.infected_agents.setGeometryLocation(location, pointMoveTo);
        }
        	
        else        	
        {
        	world.healthy_agents.setGeometryLocation(location, pointMoveTo);
        }
        if (location.geometry.getCoordinate() != c)
        {
            location.geometry.apply(pointMoveTo);     	
        }   
    }

    /** return geometry representing agent location */
    public MasonGeometry getGeometry()
    {
        return location;
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
    
    public GeomPlanarGraphEdge chooseDestination()
    // return the id of the destination edge.
    {
    	Coordinate nextBuilding = new Coordinate();

    	//determine if return home
    	Random rand_num = new Random();
    	double r = rand_num.nextDouble();
    	if (r>move_probability) {
    		//return to home node
    		nextBuilding=homeBldg;
    	} 
    	else {
    		int nextArea = (int) sampleFromDistr(area_list,area_probs);
    		ArrayList<Coordinate> bldgs_in_area = world.areasToBuildings.get(nextArea);
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
    			int num_bldgs = bldgs_in_area.size();
    			int nextPos = ThreadLocalRandom.current().nextInt(num_bldgs);
    			nextBuilding = bldgs_in_area.get(nextPos); 

    		}
    	}
    	GeomPlanarGraphEdge nextEdge = world.buildingsToEdges.get(nextBuilding);
        destination_building = nextBuilding;
        world.buildingsToVisits.put(nextBuilding, Integer.valueOf(world.buildingsToVisits.get(nextBuilding)+1));
        return nextEdge;
    }
    
    //Get a random sample from a distribution with PDF defined by xAxis and yAxis
    private double sampleFromDistr(ArrayList<Double> xAxis, ArrayList<Double> yAxis) {
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

    public Coordinate getDestinationBuilding()
    {
    	return destination_building;
    }
    
    public boolean getInBuildingStatus()
    {
    	return inBuilding;
    }
    
//    public void setMobilityAfterInfection()
//    {
//    	if (superspreader == true)
//    	{
//    		move_probability = sampleFromDistr(xAxisSuperSpreader, yAxisSuperSpreader) + 0.01; 
//    	}
//    	else
//    	{
//    		move_probability = sampleFromDistr(xAxisLowSpreader, yAxisLowSpreader) + 0.01; 
//    	}
//    }
}
