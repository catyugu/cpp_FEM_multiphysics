# Material

The `Material` class is a core component for defining the physical properties of different regions in the simulation domain.

## Overview

A `Material` object holds a collection of properties, such as thermal conductivity, electrical conductivity, or magnetic permeability. These properties can be simple constant values or complex functions that depend on other field variables (like temperature).

## Key Classes

- `Core::Material`: Represents a material and manages its properties.
- `Core::MaterialProperty`: A helper class to define properties that depend on other variables in a structured way.

## Defining Properties

There are three ways to define a material property:

### 1. Constant Values

For properties that do not change.

```cpp
auto copper = std::make_shared<Core::Material>(0, "Copper");
copper->setProperty("thermal_conductivity", 401.0); // W/(m*K)
```

### 2. Function-Based

You can provide a C++ lambda function for properties that depend on other variables.

```cpp
copper->setProperty("electrical_conductivity", [](const std::map<std::string, double>& vars) {
    double T = vars.at("Temperature");
    // Return electrical conductivity as a function of temperature T
    return 5.96e7 / (1.0 + 0.0039 * (T - 293.15));
});
```

### 3. Using `MaterialProperty` (Recommended)

The `MaterialProperty` class provides a more structured way to define variable-dependent properties.

```cpp
Core::MaterialProperty electrical_sigma(
    "electrical_conductivity",
    [](const std::map<std::string, double>& vars) -> double {
        double T = vars.at("Temperature");
        return 5.96e7 / (1.0 + 0.0039 * (T - 293.15));
    },
    {"Temperature"}  // List of dependencies
);

copper->setMaterialProperty("electrical_conductivity", electrical_sigma);
```

## Retrieving Property Values

The method for retrieving a property value depends on the context.

### In Physics Assembly (at Quadrature Points)

This is the **recommended and most accurate method**, used within the `assemble()` loop of a physics field. It ensures the highest fidelity for non-linear materials.

- `getPropertyAtQuadraturePoint(prop_name, interpolated_vars)`: Evaluates the property using variable values that have been interpolated at a specific quadrature point.

**Workflow:**
1. Inside the `assemble()` loop, at each quadrature point, use `Utils::InterpolationUtilities::interpolateAtQuadraturePoint` to get the interpolated values of required variables (e.g., "Temperature").
2. Pass these interpolated values to `material.getPropertyAtQuadraturePoint()`.

```cpp
// Inside a physics field's assemble() method
// ... looping over quadrature points ...
auto interpolated_vars = Utils::InterpolationUtilities::interpolateAtQuadraturePoint(...);
double k = material.getPropertyAtQuadraturePoint("thermal_conductivity", interpolated_vars);
```

### Element-Based Evaluation

This method evaluates properties using the variable values stored on the element (typically an average from a previous iteration). It is less accurate but can be useful for updating certain states for the element as a whole.

- `getPropertyFromElement(prop_name, element)`: Evaluates the property using the variable values associated with an element.

## API Reference

### `void setProperty(name, value)`
Sets a constant property value.

### `void setProperty(name, function)`
Sets a variable-dependent property using a lambda function.

### `void setMaterialProperty(name, material_prop)`
Sets a property using a `MaterialProperty` object.

### `double getPropertyAtQuadraturePoint(name, interpolated_vars)`
Gets the property value at a quadrature point, using interpolated field variables.

### `double getPropertyFromElement(name, element)`
Gets the property value using an element's variable values.

### `bool hasProperty(name)`
Checks if a property is defined for the material.
