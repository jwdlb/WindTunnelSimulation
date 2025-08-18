# WindTunnelSim

WindTunnelSim is a C++ simulation tool. It models a fluid within a 2D wind tunnel environemnt, around an object of the user's choice. Utilising the Cuda GPU acceleration libary, it provides a high frame rate, high resolution image.

<img src="Recources/photos/simulationImg.png" style="filter: drop-shadow(0 0 5px black);">

## Features

- **GPU acceleraion**: Utilising the **Cuda** libary to optimise simulation, through parallel proccessing.
- **Customisability**: Allows the user to choose the object in the tunnel, as well as wind speed and smoke inlet size.
- **Realistic fluid flow**: Using methods such as semi-lagrangian advection, and iterative gauss-siedel, it provides a realistic simulation for fluids that are incompressable, and non-viscous.

## Comparsion the real world

| Simulation | Aerofoil Template |
|-----------|-----------------|
| ![Simulation](Recources/photos/simulationImg4.png) | ![Aerofoil](Recources/photos/simulationImgtemplateAerofoil.jpg) |


## Prerequisites

- **Java Runtime Environment (JRE)**: Ensure you have Java installed. You can download it from the [official Oracle website](https://www.oracle.com/java/technologies/downloads/).

## Running the Application

To execute the EcoSim JAR file, download the pre-built JAR file from the [releases section](https://github.com/BRSY1/ecosim/releases):

**Using the Command Line**:
   - Open a terminal or command prompt.
   - Navigate to the directory containing `ecosim.jar`.
   - Run the following command:
     
     ```bash
     java -jar ecosim-1.0.0.jar
     ```
     
   - Ensure that the JAR file has execute permissions. If not, modify them using:
     
     ```bash
     chmod +x ecosim-1.0.0.jar
     ```
     
   - If double-clicking the JAR file doesn't launch the application, running it via the command line as shown above is recommended.

## Contributing

We welcome contributions! Please fork the repository and submit pull requests for any enhancements or bug fixes.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.






