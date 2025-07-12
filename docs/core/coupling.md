# **Core Coupling Components**

This document describes the coupling mechanism used in the C++ FEM Multiphysics project. This system is designed to be modular, allowing for different types of physics interactions to be defined and managed cleanly.

## **Classes**

### **CouplingManager**

* **Description**: Manages the coupling between physics modules. It registers all physics modules and their coupling requirements, checks if the requirements are met, and sets up the coupling.
* **Public Functions**:
  * registerField(Physics::PhysicsField& field)
  * addCoupling(std::unique\_ptr\<Coupling\> coupling)
  * setupCouplings()
* **Private Members**:
  * fields\_: std::vector\<Physics::PhysicsField\*\>
  * couplings\_: std::vector\<std::unique\_ptr\<Coupling\>\>

### **Coupling**

* **Description**: A virtual base class for all coupling types. It defines the interface for a coupling strategy.
* **Public Functions**:
  * \~Coupling(): Virtual destructor.
  * setup(std::vector\<Physics::PhysicsField\*\>& fields): Pure virtual function to set up the coupling between the registered fields.

### **ElectroThermalCoupling**

* **Description**: A concrete implementation of the Coupling class for electro-thermal coupling. It links the "Voltage" and "Temperature" fields.
* **Public Functions**:
  * setup(std::vector\<Physics::PhysicsField\*\>& fields) override