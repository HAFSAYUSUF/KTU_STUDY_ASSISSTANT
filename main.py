import os
import streamlit as st
from dotenv import load_dotenv
import google.generativeai as genai

# Load your .env file
load_dotenv()
genai.configure(api_key=os.getenv("API_KEY"))

# Initialize Gemini model
model = genai.GenerativeModel("models/gemini-1.5-flash")

# Subject-wise module content (KTU 2019 Scheme for now; 2024 to be added)
syllabus_2019 = {
    "DSP": {
        "Module 1": """Basic Elements of a DSP system, Typical DSP applications, Finite length Discrete transforms, Orthogonal transforms, 
The Discrete Fourier Transform: DFT as a linear transformation (Matrix relations), Relationship of the DFT to other transforms, IDFT, 
Properties of DFT and examples, Circular convolution, Linear Filtering methods based on the DFT (linear convolution using circular convolution), 
Filtering of long data sequences, Overlap Save and Overlap Add methods, Frequency Analysis of Signals using the DFT (concept only).""",

        "Module 2": """Efficient Computation of DFT: Fast Fourier Transform Algorithms, 
Radix-2 Decimation in Time and Decimation in Frequency FFT Algorithms, IDFT computation using Radix-2 FFT Algorithms, 
Application of FFT Algorithms – Efficient computation of DFT of Two Real Sequences and a 2N-Point Real Sequence.""",

        "Module 3": """Design of FIR Filters – Symmetric and Anti-symmetric FIR Filters, 
Design of linear phase FIR filters using Window methods (Rectangular, Hamming, Hanning), 
Frequency Sampling Method, Comparison of FIR Design Methods, 
Design of IIR Digital Filters from Analog Filters (Butterworth), IIR Filter Design using Impulse Invariance, 
IIR Filter Design using Bilinear Transformation, Frequency Transformations in the Analog and Digital Domain.""",

        "Module 4": """Structures for realization of Discrete Time Systems – Block diagram and signal flow graph representations, 
FIR Filter Structures – Direct form, cascade form, 
IIR Filter Structures – Direct form, cascade form, parallel form, 
Computational complexity of digital filters, 
Multi-rate DSP – Decimation, Interpolation, Time and Frequency domain interpretation, 
Anti-aliasing and Anti-imaging filters.""",

        "Module 5": """Computer architecture for DSP – Harvard Architecture, pipelining, MAC, 
TMS320C67xx processor – Functional block diagram, 
Finite word length effects – fixed-point vs floating-point, ADC quantization noise, 
Finite word length effects in IIR filters (coefficient quantization errors), 
Finite word length effects in FFT algorithms (round-off)."""
    },

    "LIC": {
        "Module 1": """741 Op-Amp block diagram, ideal parameters, open-loop configs, differential amplifier using BJT, current mirror, Wilson and Widlar current mirrors.""",
        "Module 2": """Negative feedback types, Op-amp feedback circuits, Inverting/Non-inverting amplifiers, Summer, Voltage follower, Instrumentation amplifier, Integrator, Differentiator, Precision rectifiers, Comparators, Schmitt trigger, Log/Antilog amplifier.""",
        "Module 3": """Op-amp based oscillators: Phase shift, Wien-bridge, waveform generators, multivibrators, active filter design: LPF, HPF, BPF, BRF, State variable filters.""",
        "Module 4": """Timer IC 555: Astable/Monostable, VCO concepts, LM566 VCO applications, PLL operation, IC 565, PLL applications.""",
        "Module 5": """Voltage regulators (fixed/adjustable), IC723 applications, protections (current limit, short circuit), DAC (R-2R, weighted), ADC (Flash, SAR)."""
    },

    "VLSI Design": {
        "Module 1": """Moore's Law, ASIC design types (Full custom, Standard cell, Gate array), SoCs, FPGA devices, Design flows (Top-Down, Bottom-Up), Logical and Physical design, Speed-power-area trade-offs.""",
        "Module 2": """Static CMOS logic design – NMOS logic (basic gate analysis), CMOS inverter (static & transient analysis), Power dissipation and delays, Static CMOS logic realization, Pass transistor logic, Transmission gate logic.""",
        "Module 3": """Dynamic logic – Precharge-Evaluate, Domino logic, NP-Domino logic, ROM – 4x4 OR/NOR/NAND, SRAM (6T CMOS), DRAM (3T/1T cells).""",
        "Module 4": """Arithmetic Circuits – Static adder, Carry-bypass, Carry-select (linear/square root), Array multiplier.""",
        "Module 5": """MOSFET Fabrication – CZ process, wafer prep, oxidation (dry/wet), diffusion/implantation, epitaxy, lithography, etching, metal deposition. Twin-tub process. Layout design – Stick diagrams, Design rules (micron/lambda), CMOS inverter/NAND/NOR layout."""
    },

    "Information Theory and Coding": {
        "Module 1": """Entropy, Properties of Entropy, Joint and Conditional Entropy, Mutual Information, Properties of Mutual Information. Discrete memoryless sources, Source code, Average length of source code, Bounds on average length, Uniquely decodable and prefix-free source codes. Kraft Inequality (with proof), Huffman code. Shannon’s source coding theorem (both achievability and converse) and operational meaning of entropy.""",
        "Module 2": """Discrete memoryless channels, Capacity of DMCs, BSC, BEC and their capacities, Channel code and rate, Shannon’s channel coding theorem and its meaning. AWGN modeling, Continuous input channels with power constraint, Differential entropy, Gaussian differential entropy, Shannon-Hartley theorem and its inferences – spectral efficiency, Shannon limit.""",
        "Module 3": """Groups, Rings, Fields, Vector spaces. Block codes and parameters, Linear block codes, Repetition and parity-check codes, Generator/parity-check matrices, Maximum likelihood decoding, Bounded distance decoding, Syndrome, Standard array.""",
        "Module 4": """Cyclic codes – Polynomial and matrix views, Systematic encoding. Overview (no decoding) of Hamming, BCH, and Reed-Solomon codes.""",
        "Module 5": """Convolutional codes – State/Trellis diagrams, Viterbi decoding. LDPC codes – Tanner graphs, Message passing for BEC."""
    },

    "Electromagnetics": {
        "Module 1": """Vector calculus: gradient, divergence, curl in Cartesian, cylindrical, spherical coordinates. Electric/magnetic fields, Coulomb’s, Gauss’s, and Ampere’s laws. Laplace/Poisson’s equations.""",
        "Module 2": """Capacitance/inductance of lines, energy stored, displacement current, continuity equation, magnetic vector potential, Maxwell’s equations, boundary conditions, wave equation.""",
        "Module 3": """Plane wave propagation in dielectrics, conductors. Attenuation, phase/group velocity, skin depth. Reflection/refraction, Snell’s law, Brewster angle.""",
        "Module 4": """Power density, Poynting theorem, polarization types. Transmission lines – parameters, equations, reflection coefficient, VSWR, input impedance.""",
        "Module 5": """Smith chart, impedance calculation, Rectangular waveguide – propagation modes, group/phase velocity, derivation & problems."""
    },

    "Disaster Management": {
        "Module 1": """Earth systems (lithosphere, atmosphere, hydrosphere, biosphere), cyclones, monsoon, disaster terms (hazard, risk, capacity, resilience), preparedness, mitigation, response.""",
        "Module 2": """Hazard/vulnerability types, assessments, disaster preparedness strategies and actions.""",
        "Module 3": """Disaster risk management and phases, risk reduction measures, prevention, response objectives/plans, international relief.""",
        "Module 4": """Stakeholder engagement, disaster communication methods/barriers, capacity building and assessment.""",
        "Module 5": """Common disasters in India, legislation, national policy, institutional frameworks, Sendai Framework (principles, actions, targets)."""
    },

    "Management": {
        "Module 1": """Introduction to management, levels of managers, classical, neo-classical and modern theories, system approach, roles of a manager.""",
        "Module 2": """Management process – planning, mission, goals, strategy, organization structures, leadership, motivation, control.""",
        "Module 3": """Productivity, decision making (certainty, risk, uncertainty), decision trees and models.""",
        "Module 4": """Project management – CPM/PERT, scheduling, crashing, redundancy, probability of completion.""",
        "Module 5": """Functional areas – HR, operations, marketing, finance, entrepreneurship, CSR, IPR, patents."""
    },

    "Control Systems": {
        "Module 1": """Control system basics – open loop and closed loop systems, feedback and its effects, types of control systems, linear/nonlinear, time-invariant/varying systems. Mathematical modeling – electrical/mechanical systems. Block diagram/signal flow graph, Mason’s gain formula, impulse response and transfer function.""",
        "Module 2": """Time domain analysis – standard test signals, time response of first and second order systems to step and ramp input, steady-state error, static error coefficients. Frequency domain specifications and correlation with time response.""",
        "Module 3": """System stability – BIBO, absolute stability, Routh-Hurwitz criterion, P/PI/PID controller effects. Root locus techniques – properties, construction, zero/pole effects.""",
        "Module 4": """Nyquist stability – fundamentals, gain/phase margin, Bode plot analysis, lag/lead compensator design using Bode plots.""",
        "Module 5": """State variable analysis – state equations, electrical/mechanical system representations, solutions using state transition matrix, controllability and observability, Kalman’s test."""
    },

    "Engineering Economics": {
        "Module 1": """Scarcity and choice, basic economic problems, PPC, firms and objectives, types of firms. Utility, law of diminishing marginal utility, demand and its determinants, law of demand, elasticity (concepts, measurement, applications), supply and law of supply, equilibrium, demand/supply shifts, consumer and producer surplus, taxation and deadweight loss.""",
        "Module 2": """Production function, law of variable proportion, economies of scale, isoquants/isocosts, producer equilibrium, expansion path, technical progress, Cobb-Douglas function, cost concepts (social, private, external, explicit, implicit, sunk), cost curves (short and long run), revenue, shutdown/break-even points.""",
        "Module 3": """Market structures – perfect/imperfect competition, monopoly (regulation), monopolistic competition, oligopoly, kinked demand, collusive oligopoly, non-price competition. Pricing – cost plus, target return, penetration, predatory, going rate, skimming.""",
        "Module 4": """Macroeconomics – circular flow, stock/flow, GDP, national income, goods classification, measurement methods. Inflation – causes, effects, control (monetary/fiscal). Business finance – bonds/shares, money/capital market, stock market, demat/trading accounts, SENSEX/NIFTY.""",
        "Module 5": """International trade – advantages/disadvantages, absolute/comparative advantage, Heckscher-Ohlin theory, balance of payments – components and concepts."""
    },
    "Antennas and Microwave Engineering": {
        "Module 1": """Antenna parameters – gain, directivity, beamwidth, effective aperture, effective height, polarization, radiation resistance and efficiency, field zones. Reciprocity, Helmholtz theorem (derivation), dipole fields/directivity/radiation resistance (short and half-wave).""",
        "Module 2": """Broadband antennas – Log periodic array (principle/design), Helical antenna types/design, Microstrip patch antenna design and feeding methods, Horn/parabolic dish antenna (E/H/Gain expressions), Inverted F mobile phone antenna.""",
        "Module 3": """Point source arrays – field of two isotropic sources, pattern multiplication, linear arrays, array factor, grating lobes, Broadside, Endfire, Dolph-Chebyshev arrays, phased arrays.""",
        "Module 4": """Microwave basics, rectangular cavity resonators (resonant freq derivation), klystrons, reflex klystrons (power/efficiency derivation), magnetrons, TWTs (slow wave, helix structure, gain derivation).""",
        "Module 5": """Microwave hybrid circuits – S-parameters, magic tees, hybrid rings, directional couplers, circulators, isolators, phase shifters. Microwave semiconductors – MESFET amplifiers, Gunn diodes (principles, modes, oscillators)."""
    },

    "Industrial Safety Engineering": {
        "Module 1": """Need for safety, productivity. Definitions: Accident, Injury, Unsafe act/condition, Dangerous occurrence, Reportable accidents. Accident causation theories. Safety organization and roles – management, supervisors, unions, govt., agencies. Safety officer duties. Safety committees.""",
        "Module 2": """PPEs – types, respiratory/non-respiratory, standards. Monitoring safety performance: frequency, severity, incidence, activity rates. Housekeeping – 5S. Work permits (hot/cold), confined space entry.""",
        "Module 3": """Construction safety – excavation, filling, underwater works, scaffolds, tunneling, blasting, demolition. NBC and IS standards. Ergonomics in construction, ergonomic hazards, MSDs and CTDs.""",
        "Module 4": """Machine safety – point of operation, guards, power presses, milling, welding safety, cutting. Material handling – techniques, equipment, wire ropes, chains, clamps maintenance.""",
        "Module 5": """Hazards – fire types, extinguishers, toxic release, Dow index, preliminary hazard analysis, HAZOP. Chemical hazards – classification, control, MSDS."""
    }
}

# Placeholder for 2024 scheme
syllabus_2024 = { 
    "physics_syllabus" :{
    "Module 1: Semiconductor Physics": """
     - Intrinsic semiconductor
- Derivation of density of electrons in conduction band and density of holes in valence band
- Intrinsic carrier concentration and its variation with temperature
- Extrinsic semiconductor (qualitative)
- Formation of p-n junction
- Fermi level in semiconductors (intrinsic and extrinsic)
- Energy band diagram of p-n junction
- Charge flow across a p-n junction (qualitative)
- Forward and reverse biased p-n junctions
- Diode equation (Derivation)
- V-I Characteristics of p-n junction
""",
    
    "Module 2: Semiconductor Devices": """
- Rectifiers: Full wave and Half wave
- Zener diode: V-I characteristics, Zener and Avalanche breakdown
- Tunnel diode: V-I characteristics
- Applications of Zener and Tunnel diodes
- Photonic devices (qualitative)
    - Photo detectors: Junction and PIN photodiodes
    - Applications
    - Solar cells: V-I characteristics, efficiency, panel stringing
    - Light Emitting Diode (LED): Applications
""",
    
    "Module 3: Superconductivity & Dielectrics": """
- Superconductivity: Transition temperature, Critical field, Meissner effect
- Type I and Type II superconductors
- Applications of superconductors
- Dielectrics:
    - Dielectric constant, Permittivity, Relative permittivity
    - Polarization: Types and relation with dielectric constant
    - Internal fields in liquids and solids
    - Clausius-Mossotti Relation
    - Dielectric loss (qualitative)
    - Dielectric breakdown (qualitative)
""",
    
    "Module 4: Laser & Fiber Optics": """
- Optical processes: Absorption, Spontaneous emission, Stimulated emission
- Properties of laser
- Principle of laser: Conditions for lasing (population inversion, pumping, metastable states)
- Laser components: Active medium, Optical resonant cavity
- Ruby laser: Construction and working
- Semiconductor laser (Qualitative)
- Applications of lasers

**Optical Fiber:**
- Principle of light propagation
- Types: Step index and Graded index fibers
- Numerical aperture (with derivation)
- Applications of optical fibers
- Fiber optic communication system (block diagram)
"""
},
"Mathematics for Electrical Science and Physical Science": {
    "Module 1": """Linear systems of equations: Gauss elimination, Row echelon form, Linear Independence: rank of a matrix. Solutions of linear systems: Existence, Uniqueness (without proof). The matrix Eigen Value Problem, Determining Eigen values and Eigen vector, Diagonalization of matrices. (Text 1: Relevant topics from sections 7.3, 7.4, 7.5, 8.1, 8.4)""",

    "Module 2": """Homogeneous linear ODEs of second order, Superposition principle, General solution. Homogeneous linear ODEs of second order with constant coefficients: Method to find general solution, solution of linear Initial Value Problem. Non-homogeneous ODEs (with constant coefficients): General solution, Particular solution by the method of undetermined coefficients (functions: ke^γx, kx^n, kcos(ωx), ksin(ωx), ke^(αx)cos(ωx), ke^(αx)sin(ωx)). Initial value Problem for Non-Homogeneous second order linear ODE. Solution by variation of parameters (Second Order). (Text 1: Sections 2.1, 2.2, 2.7, 2.10)""",

    "Module 3": """Laplace Transform, Inverse Laplace Transform, Linearity, First shifting theorem, Transform of derivatives. Solution of Initial value problems by Laplace transform (Second order linear ODE with constant coefficients, initial conditions at t=0). Unit step function, Second shifting theorem, Dirac delta function and its transform. (IVPs involving unit step function and delta excluded). Convolution theorem (without proof) and its application to inverse Laplace transforms. (Text 1: Sections 6.1 to 6.5)""",

    "Module 4": """Taylor series (without proof), Maclaurin series, Fourier series, Euler formulas, Dirichlet’s conditions, Fourier series of 2π and 2l periodic functions. Half range sine series, Half range cosine series. (Text 1: Sections 11.1, 11.2; Text 2: Section 10.8)"""
},
"Semiconductor Physics": {
    "Module 1": """Intrinsic semiconductor, Electron/hole densities, Intrinsic carrier concentration & its temperature dependence, 
    Extrinsic semiconductor (qualitative), p-n junction formation, Fermi level (intrinsic/extrinsic), Energy band diagram, 
    Forward and reverse bias, Diode equation (derivation), V-I characteristics.""",

    "Module 2": """Semiconductor Devices – Rectifiers (Half and Full wave), Zener diode (V-I, breakdown types), Tunnel diode (V-I), 
    Zener and Tunnel diode applications. Photonic devices (qualitative) – Photodetectors (Junction and PIN), Solar cells (V-I, efficiency, paneling), 
    Light Emitting Diode and its applications.""",

    "Module 3": """Superconductivity – Transition temp, Critical field, Meissner effect, Type I & II, Applications. 
    Dielectrics – Constant, Polarization, Permittivity, Relation between Polarization and Dielectric constant, 
    Types of polarization, Internal fields, Clausius-Mossotti Relation, Dielectric loss and breakdown (qualitative).""",

    "Module 4": """Laser – Absorption, Spontaneous & Stimulated emission, Laser properties, Principle, Conditions for lasing, 
    Population inversion, Pumping, Metastable states. Ruby laser (construction & working), Semiconductor laser (qualitative), Applications. 
    Fiber optics – Principle, Step/Graded index fibers, Numerical aperture (derivation), Applications, Fiber optic communication (block diagram)."""
},
"CHEMISTRY FOR INFORMATION SCIENCE AND ELECTRICALSCIENCE ": {
    "Module 1": """Electrochemical Cell – Electrode potential, Nernst equation for single electrode and full cell (with numericals), Reference electrodes – SHE and Calomel electrode (construction and working), Electrochemical series – applications. Glass electrode and pH measurement, Conductivity and digital conductivity meter, Li-ion battery and H₂-O₂ fuel cell (acid electrolyte) – construction and working.
Corrosion – Electrochemical corrosion mechanism (acidic and alkaline), Galvanic series, Corrosion control – Cathodic protection (sacrificial anode and impressed current), Electroplating of copper, Electroless plating of copper.""",

    "Module 2": """Materials for Electronic Applications – Nanomaterials classification (based on dimension and material), Synthesis methods – Sol-gel and Chemical reduction, Applications of nanomaterials.
Carbon nanostructures – Carbon Nanotubes, Fullerenes, Graphene, Carbon Quantum Dots – structure, properties, applications.
Polymers – Fire retardant polymers (Halogenated and Non-halogenated – examples only), Conducting polymers – Polyaniline and Polypyrrole (synthesis, properties, applications).
Organic electronic materials – OLED and DSSC (construction, working, applications).
Other advanced materials – materials in quantum computing, supercapacitors, and spintronics.""",

    "Module 3": """Molecular Spectroscopy and Analytical Techniques – Spectroscopy basics, types of spectra, molecular energy levels, Beer-Lambert law (with numericals).
Electronic spectroscopy – principle, types of transitions, role of conjugation, instrumentation, applications.
Vibrational spectroscopy – principle, number of vibrational modes, modes of CO₂ and H₂O, applications.
Thermal analysis – Dielectric Thermal Analysis (DETA) of polymers – working and application.
Electron Microscopy – Scanning Electron Microscope (SEM) – principle, instrumentation, and applications.""",

    "Module 4": """Environmental Chemistry – Water hardness: temporary and permanent, disadvantages, degree of hardness (numerical problems).
Water softening – ion exchange process (principle, procedure, advantages), reverse osmosis (principle and process), disinfection – chlorination (break point), ozone, UV.
Water quality parameters – DO, BOD, COD – definitions and significance.
Waste Management – Sewage treatment (primary, secondary, tertiary), trickling filter, UASB process.
E-waste – disposal methods (recycle, recovery, reuse). Climate chemistry – Greenhouse gases, ozone depletion, sustainable development and SDGs."""
},
"Engineering Graphics": {
    "Module 1": """Introduction: Importance of technical drawing in engineering, Types of lines, Dimensioning standards, BIS code of practice (No end-semester questions).
Projection of points in different quadrants. Projection of straight lines inclined to one plane and inclined to both reference planes. Trace of a line. Inclination of lines with reference planes. True length and true inclinations of lines inclined to both reference planes.""",

    "Module 2": """Projection of Simple Solids – Triangular, Rectangular, Square, Pentagonal, and Hexagonal Prisms and Pyramids, Cone and Cylinder (in simple positions).
Projection of solids including profile views. Solids with axis inclined to one reference plane and with axis inclined to both reference planes.""",

    "Module 3": """Sections of Solids – Sections of Prisms, Pyramids, Cone and Cylinder with vertical axis, cut by various planes. True shape of sections (except true shape given problems).
Development of Surfaces – Development of lateral surfaces of solids and solids cut by section planes (exclude problems with through holes).""",

    "Module 4": """Isometric Projection – Isometric scale, Isometric view and projections of Prisms, Pyramids, Cone, Cylinder, Sphere, Hemisphere and their combinations.
Computer Aided Drawing (CAD) – Introduction, Role of CAD in product design and development, Advantages. Creating 2D drawings with dimensions using suitable software. (CAD: only for internal evaluation)."""
},
"Introduction to Electrical and Electronics Engineering": {
    "Module 1": """Elementary concepts of DC electric circuits: Current and Voltage Division Rule, Relative potential. Capacitors & Inductors – V-I relations and energy stored.
Ohm's Law and Kirchhoff's laws – numerical problems. Star-delta conversion (resistive networks only) – numerical problems. 
Analysis of DC Electric circuits – Mesh current method and Node voltage method, matrix representation and solution of network equations. 
Magnetic circuits – MMF, field strength, flux density, reluctance, series and parallel magnetic circuits (qualitative).""",

    "Module 2": """Electromagnetic Induction – Faraday’s and Lenz’s law, statically and dynamically induced emf, Self and Mutual Inductance, Coefficient of Coupling (qualitative).
Alternating Current Fundamentals – Generation of AC voltage, waveforms, frequency, period, average value, RMS, form factor – numerical problems. 
AC Circuits – Phasor representation, Trig/Rect/Polar/Complex forms. Analysis of R, L, C circuits, RL, RC, RLC series circuits, Impedance, Power factor, Active/Reactive/Apparent power – numerical problems.
Three-phase systems – Generation, advantages, star and delta connections (balanced), relationships between line and phase values – numerical problems.""",

    "Module 3": """Introduction to Electronic Devices – Passive and active components.
PN junction diode – working and V-I characteristics. Zener diode – breakdown and voltage regulator. 
DC Power Supply – Block diagram, Half/Full/Bridge rectifiers, ripple factor with/without filters.
BJT – Construction, working, V-I characteristics, input/output for CE config, comparison of CE/CB/CC.
Biasing, Load line, Transistor as switch and amplifier. RC coupled amplifier – circuit and frequency response.
MOSFETs – Construction and working of N-channel and P-channel devices.""",

    "Module 4": """Modern Electronics and Applications:
Communication systems – General block diagram, Fiber optic communication system.
AM and FM – Concepts and block diagram of superheterodyne receivers.
Wired/Wireless Communication, GSM – block diagram.
Comparison of 3G, 4G, 5G, 6G technologies.
Instrumentation – Block diagram of digital multimeter, function generator, CRO and Lissajous patterns.
Applications – IoT-based smart homes, healthcare, agriculture (case studies only)."""
},
"Algorithmic Thinking with Python": {
    "Module 1": """Problem-Solving Strategies: Definition and importance, understanding multiple strategies.
Includes Trial and Error, Heuristics, Means-Ends Analysis, and Backtracking (Working backward).

The Problem-Solving Process: Computer as a model of computation. Steps: Understanding the problem,
Formulating a model, Developing an algorithm, Writing the program, Testing the program, and Evaluating the solution.

Essentials of Python Programming: Creating and using variables, numeric and string data types, math module,
Standard I/O using print and input, Python operators and precedence.""",

    "Module 2": """Algorithm and Pseudocode Representation:
Definition and purpose of pseudocode. Constructs: sequencing, selection (if-else, case), and repetition (for, while, repeat-until).
Sample problems: evaluating expressions, SI, finding max/min, grading systems, loops, factorials, largest of n numbers, etc.

Flowcharts: Symbols used – start/end, I/O, calculations, decisions, modules, loops, connectors.
Recommended tool: RAPTOR. Used only for visualizing algorithm control flow.""",

    "Module 3": """Selection and Iteration using Python:
if-else, elif, for loop, range, while loop.

Sequence Data Types: list, tuple, set, string, dictionary.
Arrays in Python using NumPy.

Decomposition and Modularisation:
Decomposing complex problems, modularisation, defining and using functions, multiple return values.

Recursion: Concept, call stack, recursion vs iteration, avoiding circularity.
Sample problems: Fibonacci, GCD, factorial, digit sum, top 3 elements (n ≥ 3).""",

    "Module 4": """Computational Approaches to Problem Solving (intro only – no analysis required):

Brute-Force: e.g., password guessing, padlocks.

Divide and Conquer: e.g., Merge Sort. Includes pros and cons.

Dynamic Programming: e.g., Fibonacci. Comparison with recursion.

Greedy Algorithm: e.g., max number of tasks in limited time. Motivation and comparison with DP.

Randomized Algorithms:
- Coupon Collector: Expected attempts to collect all coupons.
- Hat-Check Problem: Expected number of people who get their own hats back.

Concepts include motivations and characteristics of all approaches."""
},
        "Engineering Mathematics III": {
            "Module 1": "Fourier Integral, From Fourier series to Fourier Integral, Fourier Cosine and Sine integrals, Fourier Cosine and Sine Transform, Linearity, Transforms of Derivatives, Fourier Transform and its inverse, Linearity, Transforms of Derivative. (Text 1: sections 11.7, 11.8, 11.9)",
            "Module 2": "Complex Function: Limit, Continuity, Derivative, Analytic functions, Cauchy-Riemann Equations (without proof), Laplace’s Equations, Harmonic functions, Finding harmonic conjugate, Conformal mapping, Mappings of w= z^2, w=e^z, w=1/z, w=sin z. (Text 1: sections 13.3, 13.4, 17.1, 17.2, 17.4)",
            "Module 3": "Complex Integration: Line integrals in the complex plane (Definition & Basic properties), First evaluation method, Second evaluation method, Cauchy’s integral theorem (without proof) on simply connected domain, Independence of path, Cauchy integral theorem on multiply connected domain (without proof), Cauchy Integral formula (without proof). (Text 1: sections 14.1, 14.2, 14.3)",
            "Module 4": "Taylor series and Maclaurin series, Laurent series (without proof), Singularities and Zeros – Isolated Singularity, Poles, Essential Singularities, Removable singularities, Zeros of Analytic functions – Poles and Zeros, Formulas for Residues, Residue theorem (without proof), Residue Integration- Integral of Rational Functions of cosθ and sinθ. (Text 1: sections 15.4, 16.1, 16.2, 16.3, 16.4)"
        },
                "Solid State Devices": {
            "Module 1": (
                "Review of Semiconductor physics: Equilibrium and steady state conditions, "
                "Concept of effective mass and Fermi level, Density of states & Effective density of states, "
                "Equilibrium concentration of electrons and holes. Excess carriers in semiconductors: "
                "Generation and recombination mechanisms of excess carriers, quasi-Fermi levels. "
                "Carrier transport in semiconductors: Drift, conductivity and mobility, variation of mobility with temperature and doping, Hall Effect. "
                "Diffusion, Einstein relations, Poisson equations, Continuity equations, Current flow equations, Diffusion length, Gradient of quasi-Fermi level."
            ),
            "Module 2": (
                "PN junctions: Contact potential, Electrical Field, Potential and Charge distribution at the junction, Biasing and Energy band diagrams, Ideal diode equation. "
                "Bipolar junction transistor: Transistor action, Base width modulation, Current components in a BJT, Derivation of current components."
            ),
            "Module 3": (
                "Metal Semiconductor contacts: Electron affinity and work function, Ohmic and Rectifying Contacts, current voltage characteristics. "
                "Ideal MOS capacitor: band diagrams at equilibrium, accumulation, depletion and inversion, surface potential, CV characteristics, effects of real surfaces, threshold voltage, body effect. "
                "MOSFET - Drain current equation of enhancement type MOSFET (derivation) - linear and saturation region, Drain characteristics, transfer characteristics."
            ),
            "Module 4": (
                "MOSFET scaling: Need for scaling, constant voltage scaling and constant field scaling. Sub-threshold conduction in MOS. "
                "Short channel effects in MOSFETs: Channel length modulation, Drain Induced Barrier Lowering, Velocity Saturation, Threshold Voltage Variations and Hot Carrier Effects. "
                "MESFET and FinFET: Structure, operation and advantages."
            )
        },
        "Analog Circuits": {
            "Module 1": (
                "Wave Shaping Circuits: RC differentiating and integrating circuits, Analysis of First order RC low pass and high pass filter for step input - rise time, band width. "
                "Diode Clipping and clamping circuits. "
                "BJT/MOSFET Biasing: Need for biasing, DC load line, operating point, BJT biasing (CE configuration) – fixed bias & voltage divider bias (Design & analysis). MOSFET biasing."
            ),
            "Module 2": (
                "BJT Amplifiers: Design of RC coupled CE amplifier - Small signal analysis of CE amplifier using hybrid-π model (low and mid frequency). "
                "The high-frequency hybrid-π model of BJT, Miller effect, High frequency response of single stage CE amplifier, short circuit current gain, cut-off frequency fβ and unity gain bandwidth fT. "
                "MOSFET Amplifiers: Design of CS amplifier, Small signal analysis using hybrid-π model (mid frequency only), Small signal voltage gain, input & output impedance, CS stage with current source load and diode connected load. "
                "Multistage BJT Amplifiers: Types of multistage amplifiers, Effect of cascading on gain and bandwidth. Small signal voltage gain, input & output impedance of BJT cascode amplifier using hybrid-π model."
            ),
            "Module 3": (
                "Feedback amplifiers: The general feedback structure, Effect of negative feedback on gain, bandwidth, noise reduction and distortion. "
                "The four basic feedback topologies, Analysis of discrete BJT circuits in voltage-series and voltage-shunt feedback topologies - voltage gain, input and output impedance. "
                "Oscillators: Classification, criterion for oscillation, Wien bridge oscillator, Hartley and Crystal oscillator. (working principle and design equations of the circuits; analysis of Wien bridge oscillator only required)."
            ),
            "Module 4": (
                "Power amplifiers: Classification, Transformer coupled class A power amplifier, push pull class B and class AB power amplifiers, complementary-symmetry class B and Class AB power amplifiers, class C and D power amplifier - efficiency and distortion (no analysis required). "
                "Linear Voltage Regulators: Types of voltage regulators - series and shunt - working and design, load & line regulation, short circuit protection and fold back protection."
            )
        },
    
    "Logic Circuit Design": {
        "Module 1": """Introduction to digital circuits: Review of number systems representation and conversions, Arithmetic of Binary number systems, Signed and unsigned numbers, BCD.
Boolean algebra: Theorems, sum of product and product of sum - simplification, canonical forms- min term and max term, Simplification of Boolean expressions - Karnaugh map (up to 4 variables), Implementation of Boolean expressions using universal gates.""",
        
        "Module 2": """Combinational logic circuits: Half adder and Full adders, Subtractors, BCD adder, Ripple carry and carry look ahead adders, Decoders, Encoders, Code converters, Comparators, Parity generator, Multiplexers, De-multiplexers, Implementation of Boolean algebra using MUX.
Introduction to Verilog HDL – Basic language elements, Basic implementation of logic gates and combinational circuits.""",
        
        "Module 3": """Sequential Circuits: SR Latch, Flip flops - SR, JK, Master-Slave JK, D and T Flip flops. Conversion of Flip flops, Excitation table and characteristic equation. Shift registers - SIPO, SISO, PISO, PIPO and Universal shift registers. Ring and Johnsons counters. Design of Asynchronous, Synchronous and Mod N counters.""",
        
        "Module 4": """Finite state machines - Mealy and Moore models, State graphs, State assignment, State table, State reduction.
Logic Families: Electrical characteristics of logic gates (Noise margin, Fan-in, Fan-out, Propagation delay, Transition time, Power-delay product) - TTL, ECL, CMOS.
Circuit description and working of TTL and CMOS inverter, CMOS NAND and CMOS NOR gates."""
    },

    "Introduction to Artificial Intelligence and Data Science": {
        "Module 1": """Introduction to AI and Machine Learning: Basics of Machine Learning - types of Machine Learning systems - challenges in ML - Supervised learning model example - regression models - Classification model example - Logistic regression - unsupervised model example - K-means clustering. Artificial Neural Network - Perceptron - Universal Approximation Theorem (statement only) - Multi-Layer Perceptron - Deep Neural Network - demonstration of regression and classification problems using MLP.""",
        
        "Module 2": """Mathematical Foundations of AI and Data science: Role of linear algebra in Data representation and analysis – Matrix decomposition - Singular Value Decomposition (SVD) - Spectral decomposition - Dimensionality reduction technique - Principal Component Analysis (PCA).""",
        
        "Module 3": """Applied Probability and Statistics for AI and Data Science: Basics of probability - random variables and statistical measures - rules in probability - Bayes theorem and its applications - statistical estimation - Maximum Likelihood Estimator (MLE) - statistical summaries - Correlation analysis - linear correlation (direct problems only) - regression analysis - linear regression (using least square method).""",
        
        "Module 4": """Basics of Data Science: Benefits of data science - use of statistics and Machine Learning in Data Science - data science process - applications of Machine Learning in Data Science - modelling process - demonstration of ML applications in data science - Big Data and Data Science.
(For visualization the software tools like Tableau, PowerBI, R or Python can be used. For Machine Learning implementation, Python, MATLAB or R can be used.)"""
    },

    "Engineering Economics": {
        "Module 1": """Basic Economics Concepts - Basic economic problems – Production Possibility Curve – Utility – Law of diminishing marginal utility – Law of Demand - Law of supply – Elasticity - measurement of elasticity and its applications – Equilibrium - Changes in demand and supply and its effects.""",
        
        "Module 2": """Production function - Law of variable proportion – Economies of Scale – Internal and External Economies – Cobb-Douglas Production Function.""",
        
        "Module 3": """Cost concepts – Social cost, private cost – Explicit and implicit cost – Sunk cost - Opportunity cost - short run cost curves - Revenue concepts. Firms and their objectives – Types of firms – Markets - Perfect Competition – Monopoly - Monopolistic Competition - Oligopoly (features and equilibrium of a firm).""",
        
        "Module 4": """Monetary System – Money – Functions - Central Banking – Inflation - Causes and Effects – Measures to Control Inflation - Monetary and Fiscal policies – Deflation.
Taxation – Direct and Indirect taxes (merits and demerits) - GST.
National income – Concepts - Circular Flow – Methods of Estimation and Difficulties - Stock Market – Functions - Problems faced by the Indian stock market - Demat Account and Trading Account – Stock market Indicators - SENSEX and NIFTY.""",
        
        "Module 5": """Value Analysis and value Engineering - Cost Value, Exchange Value, Use Value, Esteem Value - Aims, Advantages and Application areas of Value Engineering - Value Engineering Procedure - Break-even Analysis - Cost-Benefit Analysis - Capital Budgeting - Process planning."""
    }
}

# Add 2024 subjects and modules here

# Combined scheme dictionary
scheme = {
    "2019": syllabus_2019,
    "2024": syllabus_2024
}


# Sidebar
# Streamlit UI
import streamlit as st

# Assume you already have syllabus_2019 and syllabus_2024
scheme = {
    "2019": syllabus_2019,
    "2024": syllabus_2024
}
with st.sidebar:
    st.markdown("🧠 Powered by Gemini 1.5 Flash")
    selected_year = st.selectbox("Select KTU Scheme", options=["2019", "2024"])
    st.markdown(f"📚 Syllabus: KTU ECE {selected_year} Scheme")
    st.markdown("👩‍🏫 Ask short answers, derivations, MCQs")

selected_scheme = scheme[selected_year]

# ✅ DROPDOWN FIX (ONLY FOR 2024)
if selected_scheme:
    if selected_year == "2024":
        selected_part = st.selectbox("Select Part", list(selected_scheme.keys()))
        selected_subjects = list(selected_scheme[selected_part].keys())
        subject = st.selectbox("Select Subject", selected_subjects)
        modules = list(selected_scheme[selected_part][subject].keys())
        module = st.selectbox("Select Module", modules)
        topic_content = selected_scheme[selected_part][subject][module]
    else:
        subject = st.selectbox("Select Subject", list(selected_scheme.keys()))
        modules = list(selected_scheme[subject].keys())
        module = st.selectbox("Select Module", modules)
        topic_content = selected_scheme[subject][module]

    # Show topic summary
    st.info(f"📘 {subject} - {module}: {topic_content}")

    # Initialize chat history
    if "chat_history" not in st.session_state:
        st.session_state.chat_history = []

    # Display chat history
    for role, msg in st.session_state.chat_history:
        with st.chat_message(role):
            st.markdown(msg)

    # Input from user
    user_input = st.chat_input(f"Ask from {subject} - {module}...")
    if user_input:
        st.chat_message("user").markdown(user_input)
        st.session_state.chat_history.append(("user", user_input))

        # Prompt for Gemini
        system_prompt = f"""
You are a helpful AI tutor trained in KTU ECE {selected_year} scheme.
The student is asking about **{subject} - {module}**.
Module topics: {topic_content}
Answer clearly, concisely, and in KTU exam style with formulas or short examples if needed.
"""
        full_prompt = f"{system_prompt}\n\nUser: {user_input}"
        response = model.generate_content(full_prompt)
        reply = response.text

        st.chat_message("assistant").markdown(reply)
        st.session_state.chat_history.append(("assistant", reply))
else:
    st.warning("No syllabus added for the selected scheme yet. Please update the 2024 scheme.")
