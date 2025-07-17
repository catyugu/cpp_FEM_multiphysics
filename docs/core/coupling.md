# **Core Coupling Components**

This document describes the coupling mechanism used in the C++ FEM Multiphysics project. This system is designed to be modular, allowing for different types of physics interactions to be defined and managed cleanly.

## **Classes**

### **CouplingManager**

* **Description**: Manages the coupling between physics modules. It registers all physics modules and their coupling requirements, checks if the requirements are met, and sets up the coupling.
* **Public Functions**:
  * `registerField(Physics::PhysicsField& field)`
  * `addCoupling(std::unique_ptr<Coupling> coupling)`
  * `setupCouplings()`: Sets up all registered couplings by calling their respective `setup` methods.
  * `executeCouplings()`: **New method** that iterates through all registered couplings and calls their `execute()` method, enabling bidirectional interactions between physics fields.
* **Private Members**:
  * `fields_`: `std::vector<Physics::PhysicsField*>`
  * `couplings_`: `std::vector<std::unique_ptr<Coupling>>`

### **Coupling**

* **Description**: A virtual base class for all coupling types. It defines the interface for a coupling strategy.
* **Public Functions**:
  * `~Coupling()`: Virtual destructor.
  * `setup(std::vector<Physics::PhysicsField*>& fields)`: Pure virtual function to set up the coupling between the registered fields.
  * `execute()`: **New pure virtual method** that defines the core logic for how the coupled fields interact and exchange information (e.g., calculating source terms).

### **ElectroThermalCoupling**

* **Description**: A concrete implementation of the `Coupling` class for electro-thermal coupling. It links the "Voltage" and "Temperature" fields, specifically calculating Joule heating as a volumetric source term for the heat field.
* **Public Functions**:
  * `setup(std::vector<Physics::PhysicsField*>& fields) override`: Identifies and stores pointers to the "Voltage" (electromagnetic) and "Temperature" (heat) fields from the problem's registered fields.
  * `execute() override`: Implements the coupling logic. It calculates the Joule heating (power dissipation) based on the current temperature-dependent electrical conductivity and the voltage gradient, then adds this as a `VolumetricSource` to the heat field. This method now handles dimension-specific calculations for 1D, 2D, and 3D elements.