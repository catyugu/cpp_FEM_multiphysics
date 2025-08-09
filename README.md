# C++ FEM Multiphysics Solver

`cpp-fem-multiphysics` 是一个基于现代C++开发的先进有限元分析(FEA)框架，采用模块化设计，专注于解决复杂的多物理场问题。该框架具有P级细化功能，支持电热耦合等多物理场分析，并具有高阶近似能力。

-----

## 核心特性

* **强大的有限元核心**：包含模块化的`Mesh`、`Node`、`Element`、`DOFManager`、`ElementGeometry`和`FEValues`等基础类。
* **P级细化（高阶单元）**：在几何线性网格上支持高阶数学近似（目前三角形/四面体支持到2阶，线单元支持到5阶），提高精度而无需重新划分网格。
* **高级自由度(DOF)管理**：`DOFManager`智能处理顶点和高阶（边缘）节点的自由度，支持复杂的P级细化和向量场分析。
* **多物理场支持**：
    * **热传导**：解决1D、2D和3D稳态及瞬态热传导问题。
    * **电磁场**：模拟1D、2D和3D电流分布（电压场）。
    * **磁静力学**：使用磁矢量势解决1D、2D和3D磁静力学问题。
    * **耦合物理场**：处理强耦合电热问题，其中电导率随温度变化，焦耳热作为热源。
* **现代C++**：采用C++17特性，实现简洁、高效和可维护的代码。
* **先进的求解器**：利用强大的**Eigen**库进行稀疏线性代数计算（支持`SparseLU`和`BiCGSTAB`求解器），并提供`SolverFactory`选择适当的解决策略（单场与耦合场）。
* **灵活的边界条件**：支持Dirichlet（固定值）、Neumann（通量）和Cauchy（混合/对流）边界条件，并对高阶单元进行鲁棒处理。
* **材料属性**：支持常数和**温度相关**的材料属性，使用泛化的依赖关系处理系统。
* **网格处理**：包含内置的均匀1D、2D和3D网格生成器，以及COMSOL的`.mphtxt`和**Gmsh的`.msh`**文件格式的导入器。
* **测试系统**：使用GoogleTest的综合测试套件验证各个组件和耦合模拟的正确性。
* **结果导出**：将结果导出为`.vtk`文件，以便在ParaView或VisIt等标准工具中可视化，支持节点（点）数据和不同的单元类型（线、三角形和四面体）。

-----

## 项目结构

项目组织为清晰、逻辑的目录结构：

```
cpp_FEM_multiphysics/
├── CMakeLists.txt          # 主构建脚本
├── README.md               # 当前文档
├── include/                # 头文件 (.hpp)
│   ├── core/               # 核心架构组件 (Mesh, Problem, Element)
│   │   ├── bcs/            # 边界条件类
│   │   ├── coupling/       # 耦合机制类
│   │   ├── mesh/           # 网格相关组件
│   │   └── sources/        # 源项类
│   ├── io/                 # 输入/输出工具 (VTK导出器, COMSOL导入器)
│   ├── physics/            # 特定物理场模块 (Heat2D, Current1D, Magnetic3D等)
│   ├── solver/             # 求解策略类 (SingleFieldSolver, CoupledSolver)
│   ├── post/               # 后处理工具
│   └── utils/              # 通用工具 (SimpleLogger, Quadrature, ShapeFunctions, Exceptions)
├── src/                    # 源文件 (.cpp)
│   ├── core/
│   │   ├── bcs/
│   │   ├── coupling/
│   │   ├── mesh/
│   │   └── sources/
│   ├── io/
│   ├── physics/
│   ├── solver/
│   ├── post/
│   └── utils/
├── tests/                  # 用于验证的Google Test源文件
├── docs/                   # 每个命名空间的详细markdown文档
└── data/                   # 示例网格和结果文件
```

-----

## 系统要求

* 支持C++17的编译器 (如GCC, Clang, MSVC)
* CMake >= 3.17
* Eigen >= 3.3
* GoogleTest (由CMake自动获取)

-----

## 快速开始

### 构建说明

1.  **克隆仓库:**
    ```bash
    git clone <repository-url>
    cd cpp_FEM_multiphysics
    ```

2.  **构建项目:**
    ```bash
    mkdir build && cd build
    cmake ..
    make
    ```
    这将创建两个主要的可执行文件：`cpp_FEM_multiphysics`和`run_tests`。

### 运行模拟

主可执行文件运行一个预配置的2D热传导问题模拟，使用从COMSOL文件导入的圆形网格。

运行方法：

```bash
./cpp_FEM_multiphysics
```

这将生成日志文件`femsolver.log`和结果文件`results.vtk`，您可以在ParaView中打开查看。

### 运行测试

项目包含一套完整的测试，用于验证求解器的功能。

运行测试：

```bash
./run_tests
```

所有测试应该通过，确认FEM实现对比分析解和基准测试是正确的。

-----

## 详细文档

要深入了解架构和API，请参阅每个命名空间的详细文档：

* **[Core命名空间](docs/core.md)**: 描述FEM框架的基本构建块。
* **[Core网格组件](docs/core/mesh.md)**: 详细介绍网格结构、单元以及带有`ElementGeometry`和`FEValues`的"智能单元"重构。
* **[Core耦合组件](docs/core/coupling.md)**: 解释通用耦合机制。
* **[PhysicsField命名空间](docs/physics.md)**: 详述不同物理场模拟的具体实现和P级细化策略。
* **[Solver命名空间](docs/solver.md)**: 解释可用的求解策略及其增强功能。
* **[IO命名空间](docs/io.md)**: 涵盖导入和导出数据的工具，包括更新的VTK导出功能。
* **[Utils命名空间](docs/utils.md)**: 提供通用工具如`ShapeFunctions`、`Quadrature`和`SimpleLogger`。
* **[Exceptions命名空间](docs/exceptions.md)**: 描述用于健壮错误处理的自定义异常类型。

-----

## 当前架构概述

该框架已发展成为一个复杂、模块化的系统，能够处理高阶近似的复杂多物理场模拟。当前架构围绕几个关键组件构建：

### 核心组件架构

* **Problem**: 作为中央协调器，管理物理场、耦合机制和求解过程。
* **Mesh**: 表示离散域，支持各种单元类型和细化策略。
* **Element**: 实现特定的有限元公式（三角形、四面体等），支持不同的近似阶数。
* **DOFManager**: 管理各类自由度，支持高阶近似和P级细化。
* **FEValues**: 计算并缓存形函数、其导数、雅可比行列式等在积分点的几何量。

### 物理场模块系统

框架支持各种物理场模块，每个模块都实现特定的单元公式：
* **Heat2D/Heat3D**: 用于热传导分析，支持稳态和瞬态解决方案。
* **Current2D/Current3D**: 用于电流流动分析。
* **Magnetic3D**: 用于3D磁静力学分析，采用矢量场近似。

### 耦合机制

一个灵活的系统处理多物理场相互作用：
* **Coupling**: 定义场交互接口的抽象基类。
* **ElectroThermalCoupling**: 用于耦合电热分析的实现，包含焦耳热效应。
* **SourceTerm**: 表示方程中各种物理源项的基类。

### 边界条件系统

* 支持Dirichlet、Neumann和Cauchy边界条件
* 智能合并重复约束，防止在共享边界节点上重复应用约束导致的求解失败
* 基于标签的系统，使接口更直观、更不容易出错

### 求解器和性能

* 直接稀疏求解器，混合阶耦合求解器，瞬态求解器
* 支持不同物理场具有不同单元阶数的系统稳定求解

-----

## 未来开发计划

我们正在不断扩展这个框架的功能。以下是未来开发的关键领域：

### 1. 物理场与问题定义解耦（依赖注入）

* **目标**: 消除`PhysicsField`子类对`Problem`类的直接依赖，使物理场模块独立、可测试和可重用。
* **实现策略**: 创建`MaterialManager`类或使用`Problem`类作为材料提供者，修改`PhysicsField`基类中的`assemble`方法签名以接受材料管理器的引用。

### 2. 实现耦合场分析的几何信息缓存

* **目标**: 解决非线性迭代或瞬态分析中由于重复计算不变的几何信息而导致的性能瓶颈。
* **实现策略**: 将`FEValues`转变为持久对象，实现懒加载/缓存，重构耦合更新函数。

### 3. 实现逐元迭代求解器（EBE）

* **目标**: 实现不依赖于显式组装全局稀疏矩阵的求解器，突破大规模问题的"内存墙"。
* **实现策略**: 创建新的求解器类`EBE_PCG_Solver`，实现其核心方法`matVecProduct`，使用"聚集-分散"过程。

### 4. 扩展电磁场和其他物理场模拟能力

* **目标**: 进一步扩展框架的物理建模能力。
* **实现**: 增强3D磁静力学功能，开发频域稳态求解器，进一步改进耦合瞬态求解器的稳健性和效率。

### 5. 泛化耦合场的材料属性

* **目标**: 使材料属性能够动态依赖于任何相关耦合场的解决方案。
* **实现**: 在`Material`类中开发更灵活的机制，基于一般输入场值集评估属性。
