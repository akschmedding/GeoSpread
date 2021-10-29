package sim.app.geo.covidsimulationcross;/*
 * Copyright 2011 by Mark Coletti, Keith Sullivan, Sean Luke, and
 * George Mason University Mason University Licensed under the Academic
 * Free License version 3.0
 *
 * See the file "LICENSE" for more information
 *
 * $Id$
 */

import java.awt.Color;
import javax.swing.JFrame;
import org.jfree.data.xy.XYSeries;
import sim.display.Console;
import sim.display.Controller;
import sim.display.Display2D;
import sim.display.GUIState;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.portrayal.geo.GeomPortrayal;
import sim.portrayal.geo.GeomVectorFieldPortrayal;
import sim.util.media.chart.TimeSeriesChartGenerator;

public class CovidSimulationCrossWithUI extends GUIState
{
    public Display2D display;
    public JFrame displayFrame;
    private GeomVectorFieldPortrayal roadsPortrayal = new GeomVectorFieldPortrayal(true);
    private GeomVectorFieldPortrayal tractsPortrayal = new GeomVectorFieldPortrayal(true);
    private GeomVectorFieldPortrayal agentPortrayal = new GeomVectorFieldPortrayal();
    private GeomVectorFieldPortrayal HagentPortrayal = new GeomVectorFieldPortrayal();
    TimeSeriesChartGenerator trafficChart;
    XYSeries maxSpeed;
    XYSeries avgSpeed;
    XYSeries minSpeed;

    protected CovidSimulationCrossWithUI(SimState state)
    {
        super(state);
    }

    /**
     * Main function
     * @param args
     */
    public static void main(String[] args)
    {
    	CovidSimulationCrossWithUI simple = new CovidSimulationCrossWithUI(new CovidSimulationCross(System.currentTimeMillis()));
        Console c = new Console(simple);
        c.setVisible(true);
    }

    /**
     * @return name of the simulation
     */
    public static String getName()
    {
        return "Covid Simulation";  
    }

    /**
     *  This must be included to have model tab, which allows mid-simulation
     *  modification of the coefficients
     */
    public Object getSimulationInspectedObject()
    {
        return state;
    }  // non-volatile


    /**
     * Called when starting a new run of the simulation. Sets up the portrayals
     * and chart data.
     */
    public void start()
    {
        super.start();

        CovidSimulationCross world = (CovidSimulationCross) state;

        maxSpeed = new XYSeries("Max Speed");
        avgSpeed = new XYSeries("Average Speed");
        minSpeed = new XYSeries("Min Speed");
        trafficChart.removeAllSeries();
        trafficChart.addSeries(maxSpeed, null);
        trafficChart.addSeries(avgSpeed, null);
        trafficChart.addSeries(minSpeed, null);

        state.schedule.scheduleRepeating(new Steppable()
        {
            public void step(SimState state)
            {
//                Gridlock world = (Gridlock) state;
//                double maxS = 0, minS = 10000, avgS = 0, count = 0;
//                for (Agent a : world.agentList)
//                {
//                  System.out.println("in UI agent pos = " + a.pointMoveTo.newValue);
////                    if (a.reachedDestination)
////                    {
////                        continue;
////                    }
//                    count++;
//                    double speed = Math.abs(a.speed);
//                    avgS += speed;
//                    if (speed > maxS)
//                    {
//                        maxS = speed;
//                    }
//                    if (speed < minS)
//                    {
//                        minS = speed;
//                    }
//                }
//                double time = state.schedule.time();
//                avgS /= count;
//                maxSpeed.add(time, maxS, true);
//                minSpeed.add(time, minS, true);
//                avgSpeed.add(time, avgS, true);
            }
        });

        //roadsPortrayal.setField(world.roads);
        //roadsPortrayal.setPortrayalForAll(new GeomPortrayal(Color.DARK_GRAY, 0.001, false));

        tractsPortrayal.setField(world.buildings);
        tractsPortrayal.setPortrayalForAll(new GeomPortrayal(Color.GREEN, true));

        agentPortrayal.setField(world.infected_agents);
        agentPortrayal.setPortrayalForAll(new GeomPortrayal(Color.RED, 0.00028, true));

        HagentPortrayal.setField(world.healthy_agents);
        HagentPortrayal.setPortrayalForAll(new GeomPortrayal(Color.BLUE, 0.00028, true));

        
        display.reset();
        display.setBackdrop(Color.WHITE);

        display.repaint();
    }

    /**
     * Called when first beginning a WaterWorldWithUI. Sets up the display window,
     * the JFrames, and the chart structure.
     */
    public void init(Controller c)
    {
        super.init(c);

        // make the displayer
//        display = new Display2D(1300, 600, this);
        display = new Display2D(660, 728, this);
        // turn off clipping
//        display.setClipping(false);

        displayFrame = display.createFrame();
        displayFrame.setTitle("Covid Simulation Display");
        c.registerFrame(displayFrame); // register the frame so it appears in
        // the "Display" list
        displayFrame.setVisible(true);

        display.attach(tractsPortrayal, "Buildings");
        //display.attach(roadsPortrayal, "Roads");
        display.attach(agentPortrayal, "Agents");
        display.attach(HagentPortrayal, "H Agents");

        // CHART
        trafficChart = new TimeSeriesChartGenerator();
        trafficChart.setTitle("Traffic Statistics");
        trafficChart.setYAxisLabel("Speed");
        trafficChart.setXAxisLabel("Time");
        JFrame chartFrame = trafficChart.createFrame(this);
        chartFrame.pack();
        c.registerFrame(chartFrame);
    }

    /**
     * called when quitting a simulation. Does appropriate garbage collection.
     */
    public void quit()
    {
        super.quit();

        if (displayFrame != null)
        {
            displayFrame.dispose();
        }
        displayFrame = null; // let gc
        display = null; // let gc
    }
}
