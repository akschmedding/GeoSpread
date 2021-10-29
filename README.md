# Requisite

Mason: Please put the jar files in the classpath.  https://cs.gmu.edu/~eclab/projects/mason/

**For your convenience, we provide a zip file (additional_jars.zip) with all the jar files needed to compile and run this project.**

Tips:
+ There may be other jar files missing; Please fix accordingly.
  If a jar library is missing: in contrib/geomason/pom.xml, you can find the url of the jar.
  For example:  https://mvnrepository.com/artifact/org.geoserver.extension/gs-ogr-core
  Then you can find the corresponding version of the jar on this website.

+ Directly compiling this tool using "mvn clean install" does not work because in pom.xml, some resources are pointing to the private servers in GMU (the developer of Mason and Geomason).

## Minimum Hardware Requisite

1 CPU; 30GB memory; 15GB hard disk.

# Environment and Project Set Up

+ The tested tool is complied on Windows 10, using Eclipse Java Neon IDE.

1. In Eclipse, click *File* -> *New* -> *Java Project* to create a new workspace.
1. In *New Java Project* window: enter the project name; set the root of the Eclipse workspace (project location) to *contrib/geomason*.
1. Click on *Next* to proceed to *Java Settings* page. In the *Libraries* page, click on *Add External JARs* to add all the required jar files (provided in the additional_jars.zip).
1. After creating the new workspace, right click on the project and choose "Build Path..." -> "Add External Archives (or jars)", and then choose the jars.


# Running

Open CovidSimulationAreaWithUI.java in Eclipse, and click *Run* -> *Run* to run the demo. You will see two windows pop out: 1) Covid Simulation Display, and 2) Covid Simulation Control. In Covid Simulation Control window, click the triangle button at the bottom left to start the simulation. If in Covid Simulation Display, a map is loaded and dots are moving, then the simulation is running properly.

Tips:
+ When performing the full scale of the simulation shown in the paper, we directly run the CovidSimulationArea.java without UI on the server.
+ We do not recommend running the full scale of the simulation using a laptop. Usually we only test <500 agents during development on our laptop.
+ The simulation is controlled by probability, and there is some randomness there. When running the same set of settings multiple times, we are expecting similar results but not exactly the same numbers.

## Running the simulation on a server
 1. In Eclipse (on your laptop): *File* -> *Export* -> *Runnable JAR File* -> *Next* -> *Extract required libraries into generated JAR*. You will get a runnable jar file. (Please pay attention to the java version.) Because in this way all the jars are packed into a single jar file, the generated jar file is expected to be ~1.2GB.
 1. Send the jar file to the server (using scp or putty/winscp, etc).
 1. Send the data folder *contrib/geomason/data* to the server.
 1. On the server,  please put the jar file and the data folder in the same folder, then run the jar:
  ```
  java -jar <name>.jar > <log_name>.txt
  ```


# Project Structure (Important File Locations)
+ **Scenarios** All the simulation scenarios and examples are located in *contrib/geomason/src/main/java/sim/app/geo*.
    - For each example, you can run it with or without UI.
+ **Data** The data needed to set up the parameters is in *contrib/geomason/data*.
+ **Parameters** Population and infection rate can be modified in *contrib/geomason/data/parameters.txt*. The first line corresponding to population, and the second line corresponding to infection rate. Some simulation cases also utilize a strictness of stay-at-home advisory parameter, initial infected population, an infection threshold in which mitigation measures are applied, percent vaccinated, and effectiveness of the vaccine.
+ **Distributions** Distributions are in *contrib/geomason/data/distr* and *contrib/geomason/data*. These distributions include how long patients are active after showing symptoms prior to isolation, how many patients per area, super spreaders per area, distance traveled, speed, household size, location type, and activity length. The distributions used vary depending on the simulation scenario.
+ **Maps** Maps are in *contrib/geomason/data* and *contrib/geomason/data/gangnam_seocho_new*.
