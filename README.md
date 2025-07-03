# ⚡ interEHD - ElectroHydroDynamics Two-Phase Flow Solver

[![OpenFOAM](https://img.shields.io/badge/OpenFOAM-v1812-blue.svg)](https://www.openfoam.com/)
[![License](https://img.shields.io/badge/License-GPL--3.0-red.svg)](LICENSE)
[![Status](https://img.shields.io/badge/Status-Active-green.svg)]()

> 🌊 Advanced OpenFOAM solver for electrohydrodynamics (EHD) two-phase flow simulations with iso-advection interface capturing

## 📖 Overview

**interEHDFoam** is a specialized computational fluid dynamics (CFD) solver built upon OpenFOAM's `interFoam` framework, specifically designed for simulating electrohydrodynamic (EHD) phenomena in two-phase flow systems. The solver combines:

- ⚡ **Electrohydrodynamics**: Electric field effects on fluid motion
- 🌊 **Two-Phase Flow**: Immiscible fluid interface tracking
- 🎯 **IsoAdvection**: Advanced geometric interface capturing
- 🔄 **Dynamic Meshing**: Adaptive mesh refinement capabilities

## ✨ Key Features

- 🔥 **Advanced Interface Capturing**: Uses isoAdvection method for sharp interface representation
- ⚡ **Electric Field Coupling**: Fully coupled electric field and fluid dynamics
- 🌊 **Surface Tension**: Accurate surface tension modeling
- 🔄 **Dynamic Meshing**: Support for mesh motion and topology changes
- 🎯 **PIMPLE Algorithm**: Robust pressure-velocity coupling
- 🚀 **Parallel Computing**: MPI parallelization support
- 📊 **Adaptive Refinement**: Optional dynamic mesh refinement

## 📐 Mathematical Formulation

### Governing Equations

The solver solves the following coupled system of equations:

#### 1. Continuity Equation
$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \mathbf{U}) = 0
$$

#### 2. Momentum Equation
$$
\frac{\partial (\rho \mathbf{U})}{\partial t} + \nabla \cdot (\rho \mathbf{U} \mathbf{U}) = -\nabla p + \nabla \cdot \boldsymbol{\tau} + \rho \mathbf{g} + \mathbf{F}_{\sigma} + \mathbf{F}_{E}
$$

where:
- $\mathbf{F}_{\sigma}$: Surface tension force
- $\mathbf{F}_{E}$: Electric body force

#### 3. Electric Field Equations
$$
\nabla \cdot (\sigma \nabla \phi) = 0
$$
$$
\mathbf{E} = -\nabla \phi
$$
$$
\mathbf{F}_{E} = \rho_e \mathbf{E}
$$

where:
- $\phi$: Electric potential
- $\sigma$: Electrical conductivity
- $\mathbf{E}$: Electric field
- $\rho_e$: Charge density

#### 4. Phase Fraction Transport
$$
\frac{\partial \alpha}{\partial t} + \nabla \cdot (\alpha \mathbf{U}) = 0
$$

## 🚀 Quick Start

### Prerequisites

- 🔧 **OpenFOAM v1812** or compatible version
- 💻 **GCC compiler** (recommended: GCC 7.0+)
- 🔗 **MPI library** for parallel computing
- 📊 **ParaView** for visualization (optional)

### Installation

1. **Clone the repository**
   ```bash
   cd $FOAM_RUN
   git clone <repository-url> interEHD
   cd interEHD
   ```

2. **Compile the solver**
   ```bash
   cd interEHDFoam
   wmake
   ```

3. **Verify installation**
   ```bash
   interEHDFoam -help
   ```

### 🎯 Running a Test Case

1. **Navigate to test case**
   ```bash
   cd case/
   ```

2. **Set up the case**
   ```bash
   ./Allrun.pre  # Preprocessing
   ```

3. **Run the simulation**
   ```bash
   # Serial run
   interEHDFoam
   
   # Parallel run (4 processors)
   mpirun -np 4 interEHDFoam -parallel
   ```

4. **Clean up**
   ```bash
   ./Allclean
   ```

## 📁 Project Structure

```
interEHD/
├── 📂 interEHDFoam/           # Main solver source code
│   ├── 📄 interEHDFoam.C      # Main solver file
│   ├── 📂 EHD/                # Electric field modules
│   │   ├── 📄 eleEqn.H         # Electric equation solver
│   │   ├── 📄 eleForce.H       # Electric force calculation
│   │   └── 📄 ...              # Other EHD modules
│   ├── 📄 alphaEqn.H          # Phase fraction equation
│   ├── 📄 UEqn.H              # Momentum equation
│   ├── 📄 pEqn.H              # Pressure equation
│   └── 📂 Make/               # Compilation settings
├── 📂 case/                   # Test case directory
│   ├── 📂 0/                  # Initial conditions
│   ├── 📂 constant/           # Physical properties
│   ├── 📂 system/             # Solver settings
│   └── 📄 Allrun.pre          # Preprocessing script
└── 📂 TwoPhaseFlow-master/    # Additional two-phase flow library
```

## ⚙️ Configuration

### Key Dictionary Files

#### `constant/electricProperties`
```cpp
// Electric field properties
electricModel       ohmic;
sigma1              1e-4;    // Conductivity phase 1
sigma2              1e-6;    // Conductivity phase 2
epsilonr1           80;      // Relative permittivity phase 1
epsilonr2           2;       // Relative permittivity phase 2
```

#### `system/controlDict`
```cpp
// Solver settings
application         interEHDFoam;
startFrom           startTime;
deltaT              1e-5;
maxCo               0.5;
maxAlphaCo          0.5;
```

## 🔧 Advanced Usage

### Electric Field Boundary Conditions

```cpp
// Electric potential boundary condition
phi
{
    type            fixedValue;
    value           uniform 1000;  // Voltage [V]
}

// Electric field boundary condition
E
{
    type            calculated;
    value           uniform (0 0 0);
}
```

### Adaptive Mesh Refinement

To enable adaptive mesh refinement:

```bash
interEHDFoam -overwrite
```

## 📊 Validation Cases

The solver has been validated against several benchmark cases:

- 🔵 **Droplet Deformation**: EHD-induced droplet deformation
- ⚡ **Taylor Cone**: Electrospray phenomena
- 🌊 **EHD Pumping**: Electroosmotic flow
- 💧 **Droplet Breakup**: Electric field-induced breakup

## 🤝 Contributing

We welcome contributions! Please follow these steps:

1. 🍴 Fork the repository
2. 🌿 Create a feature branch (`git checkout -b feature/amazing-feature`)
3. 💾 Commit changes (`git commit -m 'Add amazing feature'`)
4. 📤 Push to branch (`git push origin feature/amazing-feature`)
5. 🔄 Open a Pull Request

## 📚 References

1. **IsoAdvection Method**: Roenby, J., Bredmose, H., & Jasak, H. (2016). *Journal of Computational Physics*
2. **EHD Theory**: Saville, D.A. (1997). *Annual Review of Fluid Mechanics*
3. **OpenFOAM**: Weller, H.G., et al. (1998). *Computers in Physics*

## 📝 License

This project is licensed under the **GPL-3.0 License** - see the [LICENSE](LICENSE) file for details.

## 👨‍💻 Developer

**interEHD** is developed and maintained by:

- 🧑‍💻 **Developer**: seeeeeeeeeeer
- 📧 **Contact**: cotsqa@qq.com

## 🆘 Support

- 📧 **Issues**: Report bugs and request features via GitHub Issues
- 📖 **Documentation**: Check the `doc/` directory for detailed guides
- 💬 **Discussions**: Join the OpenFOAM community forums
- ✉️ **Direct Contact**: For specific questions, reach out to cotsqa@qq.com

## 🏆 Acknowledgments

- 🙏 OpenFOAM Foundation for the base framework
- 🌟 TwoPhaseFlow library contributors
- 🔬 Research community for validation cases

---

⭐ **Star this repository** if you find it useful!

🔗 **Connect with us**: [GitHub](https://github.com) | [OpenFOAM](https://www.openfoam.com/) 
