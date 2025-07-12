# **Coupling**

This document describes the coupling mechanism used in the C++ FEM Multiphysics project.

---
## **Classes**

### **CouplingManager**
* **Description**: Manages the coupling between physics modules. It registers all physics modules and their coupling requirements, checks if the requirements are met, and sets up the coupling.
* **Public Functions**:
    * `registerField(Physics::PhysicsField& field)`
    * `addCoupling(std::unique_ptr<Coupling> coupling)`
    * `setupCouplings()`
* **Private Members**:
    * `fields_`: `std::vector<Physics::PhysicsField*>`
    * `couplings_`: `std::vector<std::unique_ptr<Coupling>>`

### **Coupling**
* **Description**: A virtual base class for all coupling types.
* **Public Functions**:
    * `~Coupling()`: Virtual destructor.
    * `setup(std::vector<Physics::PhysicsField*>& fields)`: Pure virtual function to set up the coupling.

### **ElectroThermalCoupling**
* **Description**: A concrete implementation of the `Coupling` class for electro-thermal coupling.
* **Public Functions**:
    * `setup(std::vector<Physics::PhysicsField*>& fields) override`