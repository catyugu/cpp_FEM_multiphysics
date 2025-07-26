# **Core Coupling Components**

This document describes the coupling mechanism used in the C++ FEM Multiphysics project. This system is designed to be modular, allowing for different types of physics interactions to be defined and managed cleanly.

## **Classes**

### **CouplingManager**

* **Description**: Manages the coupling between physics modules. It registers all physics fields and the coupling strategies that link them.
* **Public Functions**:
  * `registerField(Physics::PhysicsField& field)`: Adds a physics field to be considered for coupling.
  * `addCoupling(std::unique_ptr<Coupling> coupling)`: Adds a specific coupling behavior (e.g., `ElectroThermalCoupling`).
  * `setupCouplings()`: Called by `Problem::setup()`, this method iterates through all registered coupling objects and calls their respective `setup` methods, allowing them to find and link the required physics fields.
  * `executeCouplings()`: Called within the solver loop, this method iterates through all registered couplings and calls their `execute()` method. This is where the actual physics interaction (e.g., calculating Joule heat) takes place.
* **Private Members**:
  * `fields_`: `std::vector<Physics::PhysicsField*>`
  * `couplings_`: `std::vector<std::unique_ptr<Coupling>>`

### **Coupling**

* **Description**: A virtual base class for all coupling types. It defines the interface for a coupling strategy.
* **Public Functions**:
  * `~Coupling()`: Virtual destructor.
  * `setup(std::vector<Physics::PhysicsField*>& fields)`: Pure virtual function to set up the coupling between the registered fields.
  * `execute()`: Pure virtual method that defines the core logic for how the coupled fields interact and exchange information. This method is responsible for calculating source terms or updating material properties based on the state of other fields.

### **ElectroThermalCoupling**

* **Description**: A concrete implementation of the `Coupling` class for electro-thermal coupling. It links the "Voltage" and "Temperature" fields.
* **Public Functions**:
  * `setup(std::vector<Physics::PhysicsField*>& fields) override`: Identifies and stores pointers to the "Voltage" (electromagnetic) and "Temperature" (heat) fields from the problem's registered fields.
  * `execute() override`: Implements the coupling logic. It iterates through every element in the mesh, calculates the Joule heating power based on the current temperature-dependent electrical conductivity and the voltage gradient, and then adds this power as a source term to the heat field's `F_coupling_` vector. This method correctly handles elements of any dimension (1D, 2D, or 3D).