import os
import streamlit as st
from dotenv import load_dotenv
import google.generativeai as genai
import fitz  # PyMuPDF for PDF rendering
import datetime

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
    "Physics for Electrical Science": {
        "Module 1: Semiconductor Physics": [
            "Intrinsic semiconductor",
            "Derivation of density of electrons in conduction band and density of holes in valence band",
            "Intrinsic carrier concentration and its variation with temperature",
            "Extrinsic semiconductor (qualitative)",
            "Formation of p-n junction",
            "Fermi level in semiconductors (intrinsic and extrinsic)",
            "Energy band diagram of p-n junction",
            "Charge flow across a p-n junction (qualitative)",
            "Forward and reverse biased p-n junctions",
            "Diode equation (Derivation)",
            "V-I Characteristics of p-n junction"
        ],
        "Module 2: Semiconductor Devices": [
            "Rectifiers: Full wave and Half wave",
            "Zener diode: V-I characteristics, Zener and Avalanche breakdown",
            "Tunnel diode: V-I characteristics",
            "Applications of Zener and Tunnel diodes",
            "Photonic devices (qualitative)",
            "Photo detectors: Junction and PIN photodiodes",
            "Solar cells: V-I characteristics, efficiency, panel stringing",
            "Light Emitting Diode (LED): Applications"
        ],
        "Module 3: Superconductivity & Dielectrics": [
            "Superconductivity: Transition temperature, Critical field, Meissner effect",
            "Type I and Type II superconductors",
            "Applications of superconductors",
            "Dielectric constant, Permittivity, Relative permittivity",
            "Polarization: Types and relation with dielectric constant",
            "Internal fields in liquids and solids",
            "Clausius-Mossotti Relation",
            "Dielectric loss (qualitative)",
            "Dielectric breakdown (qualitative)"
        ],
        "Module 4: Laser & Fiber Optics": [
            "Optical processes: Absorption, Spontaneous emission, Stimulated emission",
            "Properties and principle of laser",
            "Conditions for lasing: population inversion, pumping, metastable states",
            "Laser components: Active medium, Optical resonant cavity",
            "Ruby laser: Construction and working",
            "Semiconductor laser (Qualitative)",
            "Applications of lasers",
            "Optical Fiber: Principle, Types (Step index and Graded index)",
            "Numerical aperture (with derivation)",
            "Applications of optical fibers",
            "Fiber optic communication system (block diagram)"
        ]
    },

    "Mathematics for Electrical and Physical Science": {
        "Module 1": [
            "Gauss elimination",
            "Row echelon form",
            "Linear Independence, Rank of a matrix",
            "Solutions of linear systems: Existence, Uniqueness",
            "Eigenvalues and Eigenvectors",
            "Diagonalization of matrices"
        ],
        "Module 2": [
            "Homogeneous linear ODEs of second order",
            "Superposition principle",
            "General solution with constant coefficients",
            "Non-homogeneous ODEs: Undetermined coefficients",
            "Initial value problems (IVP)",
            "Variation of parameters"
        ],
        "Module 3": [
            "Laplace Transform and Inverse Laplace Transform",
            "Linearity, First shifting theorem",
            "Transform of derivatives",
            "Solution of IVPs using Laplace transform",
            "Unit step function, Second shifting theorem",
            "Dirac delta function and its transform",
            "Convolution theorem (application only)"
        ],
        "Module 4": [
            "Taylor and Maclaurin series",
            "Fourier series: Euler formulas, Dirichlet’s conditions",
            "Fourier series for periodic functions",
            "Half range sine and cosine series"
        ]
    },

    "Semiconductor Physics": {
        "Module 1": [
            "Intrinsic semiconductor",
            "Electron/hole densities",
            "Intrinsic carrier concentration & temperature dependence",
            "Extrinsic semiconductor (qualitative)",
            "p-n junction formation",
            "Fermi level (intrinsic/extrinsic)",
            "Energy band diagram",
            "Forward and reverse bias",
            "Diode equation (derivation)",
            "V-I characteristics"
        ],
        "Module 2": [
            "Rectifiers (Half and Full wave)",
            "Zener diode (V-I, breakdown types)",
            "Tunnel diode (V-I)",
            "Applications of Zener and Tunnel diodes",
            "Photodetectors (Junction and PIN)",
            "Solar cells (V-I, efficiency, paneling)",
            "LED and applications"
        ],
        "Module 3": [
            "Superconductivity – Transition temp, Critical field, Meissner effect",
            "Type I & II, Applications",
            "Dielectrics – Constant, Polarization, Permittivity",
            "Relation with dielectric constant",
            "Types of polarization",
            "Internal fields",
            "Clausius-Mossotti Relation",
            "Dielectric loss and breakdown (qualitative)"
        ],
        "Module 4": [
            "Laser – Absorption, Emission processes",
            "Principle, Conditions for lasing",
            "Ruby laser (construction & working)",
            "Semiconductor laser (qualitative)",
            "Applications",
            "Fiber optics – Principle, Types",
            "Numerical aperture (derivation)",
            "Applications",
            "Fiber optic communication system"
        ]
    },

    "Chemistry for Information Science and Electrical Science": {
        "Module 1": [
            "Electrochemical Cell – Nernst equation, SHE, Calomel electrode",
            "Electrochemical series – applications",
            "Glass electrode, pH measurement",
            "Digital conductivity meter",
            "Li-ion battery and H₂-O₂ fuel cell",
            "Corrosion – Mechanism, Galvanic series, Cathodic protection",
            "Electroplating of copper",
            "Electroless plating"
        ],
        "Module 2": [
            "Nanomaterials – classification, synthesis, applications",
            "Carbon nanostructures – CNT, Fullerenes, Graphene, CQDs",
            "Fire retardant polymers",
            "Conducting polymers – Polyaniline, Polypyrrole",
            "Organic electronics – OLED, DSSC",
            "Quantum computing, Supercapacitors, Spintronics"
        ],
        "Module 3": [
            "Spectroscopy – Basics, Beer-Lambert Law",
            "Electronic spectroscopy – principle, instrumentation",
            "Vibrational spectroscopy – modes of CO₂, H₂O",
            "DETA for polymers",
            "Electron Microscopy – SEM"
        ],
        "Module 4": [
            "Water hardness – types, numericals",
            "Water softening – ion exchange, reverse osmosis",
            "Disinfection methods",
            "Water quality parameters – DO, BOD, COD",
            "Sewage treatment – primary to tertiary",
            "E-waste disposal – recycle, reuse",
            "Climate chemistry – Greenhouse gases, SDGs"
        ]
    },

    "Engineering Graphics": {
        "Module 1": [
            "Importance of drawing",
            "Types of lines, dimensioning, BIS code",
            "Projection of points and straight lines",
            "Inclination with reference planes, True length"
        ],
        "Module 2": [
            "Projection of solids – prisms, pyramids, cone, cylinder",
            "Profile views, Axis inclined to one or both planes"
        ],
        "Module 3": [
            "Sections of solids by various planes",
            "True shape of sections",
            "Development of surfaces"
        ],
        "Module 4": [
            "Isometric projection and scale",
            "Projections of solids and combinations",
            "CAD – basics, 2D drawings with dimensions"
        ]
    },

    "Introduction to Electrical and Electronics Engineering": {
        "Module 1": [
            "DC circuits – Current/Voltage division, Relative potential",
            "V-I relations, Energy in Capacitors/Inductors",
            "Ohm's & Kirchhoff's Laws",
            "Star-delta conversion",
            "Mesh and Node methods",
            "Magnetic circuits"
        ],
        "Module 2": [
            "Faraday’s & Lenz’s laws",
            "Self and Mutual Inductance",
            "AC – Generation, RMS, Average, Form factor",
            "Phasors, RLC circuits, Power factor",
            "Three-phase systems"
        ],
        "Module 3": [
            "Passive and active components",
            "Diodes – PN and Zener",
            "Power supplies – rectifiers, ripple factor",
            "BJT – CE/CB/CC, Amplifiers, Load line",
            "MOSFETs – N and P-channel"
        ],
        "Module 4": [
            "Communication – Block diagrams",
            "AM/FM and superheterodyne receivers",
            "GSM, 3G to 6G comparison",
            "Instruments – DMM, CRO, Function generator",
            "IoT case studies – Smart homes, health, agri"
        ]
    },

    "Algorithmic Thinking with Python": {
        "Module 1": [
            "Problem-solving strategies – Trial and Error, Heuristics, etc.",
            "Problem-solving process: Understand, Model, Algorithm, Code, Test",
            "Python basics – variables, datatypes, math module, I/O, operators"
        ],
        "Module 2": [
            "Pseudocode: sequencing, selection, iteration",
            "Sample problems",
            "Flowcharts – symbols and RAPTOR tool"
        ],
        "Module 3": [
            "Selection and iteration in Python",
            "Sequence types: list, tuple, set, dict",
            "Arrays using NumPy",
            "Functions and modularisation",
            "Recursion – concept and examples"
        ],
        "Module 4": [
            "Brute-Force approach",
            "Divide and Conquer – Merge Sort",
            "Dynamic Programming – Fibonacci",
            "Greedy Algorithm",
            "Randomized Algorithms – Coupon Collector, Hat-Check Problem"
        ]
    },

    
    "Engineering Mathematics III": {
        "Module 1": ["Fourier Integral, From Fourier series to Fourier Integral, Fourier Cosine and Sine integrals, Fourier Cosine and Sine Transform, Linearity, Transforms of Derivatives, Fourier Transform and its inverse, Linearity, Transforms of Derivative. (Text 1: sections 11.7, 11.8, 11.9)"],
        "Module 2": ["Complex Function: Limit, Continuity, Derivative, Analytic functions, Cauchy-Riemann Equations (without proof), Laplace’s Equations, Harmonic functions, Finding harmonic conjugate, Conformal mapping, Mappings of w= z^2, w=e^z, w=1/z, w=sin z. (Text 1: sections 13.3, 13.4, 17.1, 17.2, 17.4)"],
        "Module 3":["Complex Integration: Line integrals in the complex plane (Definition & Basic properties), First evaluation method, Second evaluation method, Cauchy’s integral theorem (without proof) on simply connected domain, Independence of path, Cauchy integral theorem on multiply connected domain (without proof), Cauchy Integral formula (without proof). (Text 1: sections 14.1, 14.2, 14.3)"],
        "Module 4": ["Taylor series and Maclaurin series, Laurent series (without proof), Singularities and Zeros – Isolated Singularity, Poles, Essential Singularities, Removable singularities, Zeros of Analytic functions – Poles and Zeros, Formulas for Residues, Residue theorem (without proof), Residue Integration- Integral of Rational Functions of cosθ and sinθ. (Text 1: sections 15.4, 16.1, 16.2, 16.3, 16.4)"]
        },
        "Solid State Devices": {
            "Module 1": [
                "Review of Semiconductor physics: Equilibrium and steady state conditions, "
                "Concept of effective mass and Fermi level, Density of states & Effective density of states, "
                "Equilibrium concentration of electrons and holes. Excess carriers in semiconductors: "
                "Generation and recombination mechanisms of excess carriers, quasi-Fermi levels. "
                "Carrier transport in semiconductors: Drift, conductivity and mobility, variation of mobility with temperature and doping, Hall Effect. "
                "Diffusion, Einstein relations, Poisson equations, Continuity equations, Current flow equations, Diffusion length, Gradient of quasi-Fermi level."
            ],
            "Module 2": [
                "PN junctions: Contact potential, Electrical Field, Potential and Charge distribution at the junction, Biasing and Energy band diagrams, Ideal diode equation. "
                "Bipolar junction transistor: Transistor action, Base width modulation, Current components in a BJT, Derivation of current components."
            ],
            "Module 3": [
                "Metal Semiconductor contacts: Electron affinity and work function, Ohmic and Rectifying Contacts, current voltage characteristics. "
                "Ideal MOS capacitor: band diagrams at equilibrium, accumulation, depletion and inversion, surface potential, CV characteristics, effects of real surfaces, threshold voltage, body effect. "
                "MOSFET - Drain current equation of enhancement type MOSFET (derivation) - linear and saturation region, Drain characteristics, transfer characteristics."
            ],
            "Module 4": [
                "MOSFET scaling: Need for scaling, constant voltage scaling and constant field scaling. Sub-threshold conduction in MOS. "
                "Short channel effects in MOSFETs: Channel length modulation, Drain Induced Barrier Lowering, Velocity Saturation, Threshold Voltage Variations and Hot Carrier Effects. "
                "MESFET and FinFET: Structure, operation and advantages."
            ],
        },
        "Analog Circuits": {
            "Module 1": [
                "Wave Shaping Circuits: RC differentiating and integrating circuits, Analysis of First order RC low pass and high pass filter for step input - rise time, band width. "
                "Diode Clipping and clamping circuits. "
                "BJT/MOSFET Biasing: Need for biasing, DC load line, operating point, BJT biasing (CE configuration) – fixed bias & voltage divider bias (Design & analysis). MOSFET biasing."
            ],
            "Module 2": [
                "BJT Amplifiers: Design of RC coupled CE amplifier - Small signal analysis of CE amplifier using hybrid-π model (low and mid frequency). "
                "The high-frequency hybrid-π model of BJT, Miller effect, High frequency response of single stage CE amplifier, short circuit current gain, cut-off frequency fβ and unity gain bandwidth fT. "
                "MOSFET Amplifiers: Design of CS amplifier, Small signal analysis using hybrid-π model (mid frequency only), Small signal voltage gain, input & output impedance, CS stage with current source load and diode connected load. "
                "Multistage BJT Amplifiers: Types of multistage amplifiers, Effect of cascading on gain and bandwidth. Small signal voltage gain, input & output impedance of BJT cascode amplifier using hybrid-π model."
            ],
            "Module 3": [
                "Feedback amplifiers: The general feedback structure, Effect of negative feedback on gain, bandwidth, noise reduction and distortion. "
                "The four basic feedback topologies, Analysis of discrete BJT circuits in voltage-series and voltage-shunt feedback topologies - voltage gain, input and output impedance. "
                "Oscillators: Classification, criterion for oscillation, Wien bridge oscillator, Hartley and Crystal oscillator. (working principle and design equations of the circuits; analysis of Wien bridge oscillator only required)."
            ],
            "Module 4": [
                "Power amplifiers: Classification, Transformer coupled class A power amplifier, push pull class B and class AB power amplifiers, complementary-symmetry class B and Class AB power amplifiers, class C and D power amplifier - efficiency and distortion (no analysis required). "
                "Linear Voltage Regulators: Types of voltage regulators - series and shunt - working and design, load & line regulation, short circuit protection and fold back protection."
            ]
        },
    
        "Logic Circuit Design": {
            "Module 1": ["Introduction to digital circuits: Review of number systems representation and conversions, Arithmetic of Binary number systems, Signed and unsigned numbers, BCD."
"Boolean algebra: Theorems, sum of product and product of sum - simplification, canonical forms- min term and max term, Simplification of Boolean expressions - Karnaugh map (up to 4 variables), Implementation of Boolean expressions using universal gates."],
        
            "Module 2": ["Combinational logic circuits: Half adder and Full adders, Subtractors, BCD adder, Ripple carry and carry look ahead adders, Decoders, Encoders, Code converters, Comparators, Parity generator, Multiplexers, De-multiplexers, Implementation of Boolean algebra using MUX."
"Introduction to Verilog HDL – Basic language elements, Basic implementation of logic gates and combinational circuits."],
        
            "Module 3": ["""Sequential Circuits: SR Latch, Flip flops - SR, JK, Master-Slave JK, D and T Flip flops. Conversion of Flip flops, Excitation table and characteristic equation. Shift registers - SIPO, SISO, PISO, PIPO and Universal shift registers. Ring and Johnsons counters. Design of Asynchronous, Synchronous and Mod N counters."""],
        
            "Module 4": ["""Finite state machines - Mealy and Moore models, State graphs, State assignment, State table, State reduction.
Logic Families: Electrical characteristics of logic gates (Noise margin, Fan-in, Fan-out, Propagation delay, Transition time, Power-delay product) - TTL, ECL, CMOS.
Circuit description and working of TTL and CMOS inverter, CMOS NAND and CMOS NOR gates."""]
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
    },
       # Removed S3 syllabus content as requested



"Mathematics for Electrical Science": {
    "Module 1": "Random variables, Discrete random variables and their probability distributions, Cumulative distribution function, Expectation, Mean and variance, Binomial distribution, Poisson distribution, Poisson distribution as a limit of the binomial distribution, Joint pmf of two discrete random variables, Marginal pmf, Independent random variables, Expected value of a function of two discrete variables.",
    "Module 2": "Continuous random variables and their probability distributions, Cumulative distribution function, Expectation, Mean and variance, Uniform, Normal and Exponential distributions, Joint pdf of two Continuous random variables, Marginal pdf, Independent random variables, Expectation value of a function of two continuous variables.",
    "Module 3": "Confidence Intervals, Confidence Level, Confidence Intervals and One-side confidence intervals for a Population Mean for large and small samples (normal distribution and t-distribution), Hypotheses and Test Procedures, Type I and Type II error, z Tests for Hypotheses about a Population Mean (for large sample), t Test for Hypotheses about a Population Mean (for small sample), Tests concerning a population proportion for large and small samples.",
    "Module 4": "Random process concept, classification of process, Methods of Description of Random process, Special classes, Average Values of Random Process, Stationarity- SSS, WSS, Autocorrelation functions and its properties, Ergodicity, Mean-Ergodic Process, Mean-Ergodic Theorem, Correlation Ergodic Process, Distribution Ergodic Process."
},

"Signals and Systems": {
    "Module 1": "Introduction to signals and systems: Continuous time and discrete time signals - Elementary signals, Classification of signals, Basic signal operations. Continuous time and discrete time systems – Representation and Classification (memory, causal, stable, linear, time-invariant, invertible). Convolution integral and convolution sum operations. Continuous time and discrete time LTI systems - Stability and causality of LTI systems.",
    "Module 2": "Frequency domain representation of continuous time signals: Continuous time Fourier series - Exponential Fourier series representation of periodic signals. Continuous time Fourier transform - Convergence and Gibbs phenomenon, Continuous time Fourier transform of standard signals, Properties of Continuous time Fourier transform, Inverse Transform. Bilateral Laplace Transform, Concept of ROC, Relation of Laplace transform to Fourier Transform.",
    "Module 3": "Sampling of continuous time signals to discrete signals and frequency domain representation of discrete time signals: Conversion of continuous time signal to discrete time signal, Sampling theorem for low pass signals, Nyquist criteria, Aliasing. Discrete time Fourier series for discrete periodic signals. Discrete time Fourier transform (DTFT) - Convergence condition, DTFT of standard signals, Properties of DTFT, Inverse transform. Z transform - ROC, Properties (Proof not needed), Inverse transform, Relation between DTFT and Z-Transform.",
    "Module 4": "Analysis of LTI systems using Transforms: Concept of transfer function - Frequency response, Magnitude response and phase response. Analysis of Continuous time LTI systems using Laplace and Fourier transforms. Analysis of discrete time LTI systems using DTFT and Z transforms, Stability and causality using Z transform."
},

"Linear Integrated Circuits": {
    "Module 1": "Differential Amplifiers: Differential amplifier configurations using BJT, DC Analysis - transfer characteristics; AC analysis - differential and common mode gains, CMRR, input and output resistance, voltage gain, constant current bias, constant current source. Concept of current mirror: two-transistor current mirror, Wilson and Widlar current mirrors. Operational amplifiers (Op Amps): The 741 Op Amp, Block diagram, Ideal Op Amp parameters, typical parameter values for 741, equivalent circuit, open loop configurations, voltage transfer curve, frequency response curve.",
    "Module 2": "Op Amp with negative feedback: General concept of Voltage Series, Voltage Shunt, Current Series and Current Shunt negative feedback, Op Amp circuits with Voltage Series and Voltage Shunt feedback, Virtual ground concept. Analysis of inverting and non-inverting amplifier for closed loop gain, Input Resistance and Output Resistance.",
    "Module 3": "Op Amp applications: Summer, Voltage Follower, Differential and Instrumentation Amplifiers, Voltage to Current and Current to Voltage converters, Integrator, Differentiator, Precision Rectifiers, Comparators, Schmitt Triggers, Log and Antilog amplifiers.",
    "Module 4": "Differential Amplifiers (Repeat for emphasis): Differential amplifier configurations using BJT, DC Analysis - transfer characteristics; AC analysis - differential and common mode gains, CMRR, input and output resistance, voltage gain, constant current bias, constant current source. Concept of current mirror: two-transistor current mirror, Wilson and Widlar current mirrors. Operational amplifiers (Op Amps): The 741 Op Amp, Block diagram, Ideal Op Amp parameters, typical parameter values for 741, equivalent circuit, open loop configurations, voltage transfer curve, frequency response curve."

# Add more subjects here as needed

},
"Mathematics for Electrical Science": {
    "Module 1": "Random variables, Discrete random variables and their probability distributions, Cumulative distribution function, Expectation, Mean and variance, Binomial distribution, Poisson distribution, Poisson distribution as a limit of the binomial distribution, Joint pmf of two discrete random variables, Marginal pmf, Independent random variables, Expected value of a function of two discrete variables.",
    "Module 2": "Continuous random variables and their probability distributions, Cumulative distribution function, Expectation, Mean and variance, Uniform, Normal and Exponential distributions, Joint pdf of two Continuous random variables, Marginal pdf, Independent random variables, Expectation value of a function of two continuous variables.",
    "Module 3": "Confidence Intervals, Confidence Level, Confidence Intervals and One-side confidence intervals for a Population Mean for large and small samples (normal distribution and t-distribution), Hypotheses and Test Procedures, Type I and Type II error, z Tests for Hypotheses about a Population Mean (for large sample), t Test for Hypotheses about a Population Mean (for small sample), Tests concerning a population proportion for large and small samples.",
    "Module 4": "Random process concept, classification of process, Methods of Description of Random process, Special classes, Average Values of Random Process, Stationarity- SSS, WSS, Autocorrelation functions and its properties, Ergodicity, Mean-Ergodic Process, Mean-Ergodic Theorem, Correlation Ergodic Process, Distribution Ergodic Process."
},

"Microcontrollers": {
    "Module 1": "Microcontroller Architecture – General internal architecture, Address bus, Data bus, control bus. The Microcontroller 8051: Features, Block diagram, PSW, accumulator, program counter, Memory organization – RAM & ROM, register banks and stack, SFRs, I/O port organization, Interrupts.",
    "Module 2": "Instruction Set of 8051 & Addressing modes: Classification - Data transfer, arithmetic, logical, branching groups. Addressing modes - Types, Accessing data from internal and external memory.",
    "Module 3": "Programming 8051 Using Assembly: Data types & directives, Subroutine concept, Software delay programming. Programming in Embedded C: Introduction and advantages.",
    "Module 4": "Timers and Counters: Timer0, Timer1, configuration and mode programming. Serial Communication: RS-232, serial port programming – transmitting and receiving. Interrupts: External, timer and serial port interrupts, priority settings."
},

"Measurement and Instrumentation": {
    "Module 1": "Introduction to instruments, functional elements of instrumentation systems, measurement systems need and classification, Static and Dynamic characteristics of instruments.",
    "Module 2": "Sensors and Transducers – classification and selection criteria, operation of resistive (potentiometers, strain gauges, thermistors), inductive (LVDT), capacitive transducers (capacitive microphone).",
    "Module 3": "Electronic Measuring Instruments – DSO, waveform analyser, digital frequency meter, harmonic distortion meter, spectrum analyser, logic state analyser, GPIB. EMI, grounding and shielding.",
    "Module 4": "PLC Programming – Basics, Ladder diagrams, Timers and Counters, Arithmetic and data handling functions."
},

"Power Electronics": {
    "Module 1": "Introduction, Properties of ideal switch, BJT, MOSFET, IGBT – structure, characteristics and comparison, SiC and GaN basics, SOA and drive circuits.",
    "Module 2": "SCR structure, characteristics, Rectifiers – 1-phase half controlled (R load), fully controlled bridge, 3-phase half wave controlled, output equations and problems.",
    "Module 3": "DC-DC converters: Buck, Boost, Buck-boost – waveforms, voltage/current ripple, Isolated converters: Flyback, Forward, Push Pull, Half bridge, Full bridge – principles and equations.",
    "Module 4": "DC-AC Inverters – Push-Pull, Half bridge, Full bridge, PWM inverters (single and sinusoidal pulse width), RMS output equation and waveforms."
},

"Machine Learning": {
    "Module 1": "Review of supervised, unsupervised learning, PCA, SVD, instance-based vs model-based, hyperparameters, regularization, training types, data issues, overfitting/underfitting, bias/variance, performance metrics.",
    "Module 2": "Regression – linear, logistic, error functions (MSE, L1, L2, Cross entropy), multivariate regression. Classification – Naive Bayes, SVM, Decision trees, random forests, ensemble methods (boosting, bagging).",
    "Module 3": "Unsupervised learning – K-means, hierarchical clustering, distance metrics (Euclidean, Manhattan, Minkowski, cosine). Reinforcement learning – Q-learning, HMMs.",
    "Module 4": "Neural networks – biological neuron, perceptron, XOR problem, MLP, gradient descent, activation functions (Sigmoid, ReLU, tanh), backpropagation, L1/L2 regularization."
},

"Object Oriented Programming": {
    "Module 1": "Introduction to OOP, Java basics – buzzwords, structure, compiler, JVM, lexical issues.",
    "Module 2": "Java Fundamentals – primitive types, arrays, strings, operators, control statements.",
    "Module 3": "OOP in Java – classes, objects, methods, constructors, inheritance, method overloading/overriding.",
    "Module 4": "Packages, exception handling, I/O (console, file), Swing GUI, JDBC."
},

"Digital System Design": {
    "Module 1": "CSSN – analysis, modeling, state reduction and assignment, design.",
    "Module 2": "ASM charts and realization. ASC – analysis, flow table reduction, races, transition table.",
    "Module 3": "Hazards – static/dynamic, essential hazards, synchronizers, practical issues (skew, jitter), fault detection (fault table, path sensitization, Boolean difference).",
    "Module 4": "VLSI Design Flow – schematic, VHDL modeling (dataflow, behavioral, structural), synthesis and simulation."
},

"Digital Systems and VLSI Design": {
    "Module 1": "CSSN, Mealy/Moore models, modeling, state assignment and reduction, ASM charts.",
    "Module 2": "ASC – analysis, design, ALU design.",
    "Module 3": "Hazards, synchronizers, flip-flops, clock issues, fault modeling (Kohavi, BIST).",
    "Module 4": "VLSI design – FSM, HDL, VHDL processes, modeling, simulation, synthesis."
},

"Economics for Engineers": {
    "Module 1": "Basic economics concepts – demand/supply, elasticity, production function, utility, economies of scale.",
    "Module 2": "Costs – social/private, explicit/implicit, sunk/opportunity, short-run curves. Firms, market types – perfect, monopoly, oligopoly.",
    "Module 3": "Money – functions, central banking, inflation/deflation, monetary/fiscal policy, taxation (direct/indirect), GST, national income.",
    "Module 4": "Stock market – demat/trading accounts, indicators (SENSEX, NIFTY). Value Engineering – value types, procedure, breakeven, cost-benefit, capital budgeting."
},
"Engineering Ethics and Sustainable Development": {
    "Module 1": "Fundamentals of ethics - Personal vs. professional ethics, Civic virtue, Respect for others, Profession and Professionalism, Integrity in design, development, and research, Plagiarism, legal challenges, case studies. Cyberethics - data, information, and knowledge, Cybertrust, cybersecurity, data management, Social impacts of high-tech, Managing conflict, Collective bargaining, Confidentiality, Codes of Ethics.",
    
    "Module 2": "Gender Studies – basic concepts: sex, gender, sexuality, identity, gender stereotypes. Gender disparity in education, work and life, History of women in S&T, Gendered technologies, Ethical practices: equity, diversity, gender justice, Policies and empowerment for women and transgender people.",
    
    "Module 3": "Environmental Ethics – definition, development, and theories (anthropocentrism, biocentrism, ecocentrism). Sustainable Engineering – triple bottom line, life cycle analysis, sustainability metrics. Ecosystems and biodiversity – functions, loss, conservation, Indian case studies. Landscape and Urban Ecology – impacts of urbanization, sustainable urban planning.",
    
    "Module 4": "Hydrology & Water – water cycle, scarcity, pollution, sustainable practices, environmental flows. Zero Waste & Circular Economy – waste reduction, reuse/recycling, degrowth strategies. Mobility – transportation impacts, sustainable urban mobility, E-Mobility. Renewable Energy – solar, wind, hydro, biomass. Climate change – science, effects, engineering mitigation. Environmental policies – laws, compliance, ethical engineering. Case studies and future directions."
},
"Electromagnetics": {
    "Module 1": "Coordinate systems (rectangular, cylindrical, spherical), vector calculus (curl, divergence, gradient), Coulomb’s law, Gauss’s law, Ampere’s law, capacitance and inductance of two-wire and coaxial cables, magnetic scalar and vector potentials, Poisson and Laplace equations.",
    "Module 2": "Maxwell’s equations from fundamental laws, boundary conditions, wave equations and their solutions, plane wave propagation in perfect dielectrics, lossy media, good conductors, skin depth, polarization of waves.",
    "Module 3": "Reflection and refraction of EM waves (normal & oblique incidence), Snell’s law, Brewster angle, power density, Poynting vector theorem.",
    "Module 4": "Transmission line circuit model (L and C), transmission line equations, characteristic impedance, reflection coefficient, VSWR, input impedance, Smith chart calculations, rectangular waveguides – TE and TM modes, dominant mode, group and phase velocities."
},

"Analog and Digital Communication": {
    "Module 1": "Communication system block diagram, need for modulation, AM (spectrum and equation), DSB-SC, SSB, VSB, angle modulation (FM and PM), FM spectra, Carson’s rule, pre-emphasis/de-emphasis, FM receiver block diagram, superheterodyne receiver, noise types.",
    "Module 2": "Sampling, quantization, SQNR, companding, PCM (Tx & Rx), DPCM, delta modulation, slope overload, line codes.",
    "Module 3": "Baseband transmission over AWGN, ISI modeling and Nyquist criterion, raised cosine filter, equalization (zero forcing), signal space and Gram-Schmitt orthogonalization, matched filter, MAP and ML receivers.",
    "Module 4": "Bandpass modulation schemes – BPSK, QPSK, QAM – signal constellation, transmitter & receiver block diagrams, BER expressions, BER vs SNR plots."
},

"Control Systems": {
    "Module 1": "Basics of control systems – open and closed loop, modeling of electrical and mechanical translational systems, transfer function, block diagram reduction, signal flow graph and Mason’s gain formula.",
    "Module 2": "Time domain analysis – response of 1st and 2nd order systems, damping cases, time domain specifications, steady-state error, static error coefficients.",
    "Module 3": "Stability – BIBO, absolute stability, Routh-Hurwitz, root locus – properties and construction, frequency domain analysis – gain/phase margins, Bode and Nyquist plots, controllers – P, PI, PID, compensators – lag and lead (design steps only).",
    "Module 4": "State space modeling, transfer function from state equations, solution using STM, controllability and observability – Kalman’s test."
},

"Digital Signal Processing": {
    "Module 1": "Sampling, Z-transform, DTFT, DFT as matrix transformation, IDFT, DFT properties, circular vs linear convolution, overlap-add and overlap-save methods, frequency analysis using DFT.",
    "Module 2": "FIR filters – symmetric/antisymmetric, linear phase using windows (rectangular, Hamming, Hanning), IIR filters – analog prototypes (Butterworth), impulse invariance, bilinear transform, frequency transformations.",
    "Module 3": "Filter realizations – FIR (direct), IIR (direct, transposed, cascade, parallel), multirate DSP – decimation, interpolation, anti-aliasing and anti-imaging filters.",
    "Module 4": "Efficient DFT computation using FFT (Radix-2 DIT), DSP architecture – Harvard, pipelining, MAC, TMS320C67xx overview, fixed/floating point, ADC quantization noise."
},

"Biomedical Engineering": {
    "Module 1": "EMG – electrical activity of muscles, signal acquisition, myoelectric control, electrical stimulation, applications.",
    "Module 2": "Clinical lab instruments – oxymeters, cell counter, flame photometer, spectrophotometer; therapeutic devices – pacemakers, defibrillators, heart-lung machine, dialysis, surgical diathermy, ventilators.",
    "Module 3": "Biomedical telemetry – single channel telemetry system for ECG, system components and applications.",
    "Module 4": "Medical imaging – X-ray, CT (principles, reconstruction), ultrasound (A/B/M scans), real-time imaging, MRI – NMR components and biological effects."
},
    "Data Structures": {
        "Module 1": "Basic Concepts of Data Structures: Algorithms, Performance Analysis, Space Complexity, Time Complexity, Asymptotic Notations. Arrays: Linear Search and Binary Search, Stacks, Queues-Circular Queues, Priority Queues, Double Ended Queues, Evaluation of Expressions.",
        "Module 2": "Linked List: Self-Referential Structures, Dynamic Memory Allocation, Singly Linked List- Operations on Linked List. Doubly Linked List, Circular Linked List, Stacks and Queues using Linked List, Polynomial representation using Linked List.",
        "Module 3": "Trees and Graphs: Trees, Binary Trees-Tree Operations, Binary Tree Representation, Tree Traversals, Binary Search Trees- Binary Search Tree Operations. Graphs, Representation of Graphs, Depth First Search and Breadth First Search on Graphs, Applications of Graphs.",
        "Module 4": "Sorting and Hashing: Sorting Techniques – Selection Sort, Insertion Sort, Quick Sort, Merge Sort and Heap Sort. Hashing- Hashing Techniques, Collision Resolution, Overflow handling, Hashing functions – Mid square, Division, Folding, Digit Analysis."
    },

    "Sensors and Actuators": {
        "Module 1": "Introduction to Sensors and actuators, Block diagram of a closed loop control System, Sensors and Transducers, Sensor Classification, Characteristics: Transfer Function, Calibration, Span, FSO, Accuracy, Precision, Hysteresis, Nonlinearity, Saturation, Repeatability, Dead Band, Sensitivity, Resolution.",
        "Module 2": "Position and Displacement Sensors - Potentiometric, Capacitive, LVDT, Hall Effect. Pressure Sensors - Mercury, Bellows, Membranes, Thin plates, Piezoresistive, Capacitive. Force/Strain/Tactile Sensors - Strain Gauges, Piezoelectric, Capacitive Touch, Acoustic, Optical Sensors.",
        "Module 3": "Flow Sensors - Ultrasonic, Electromagnetic. Temperature Sensors - RTDs, Thermistors, Thermocouples. Proximity Sensors - PIR, Ultrasonic. Smart Sensors - Block Diagram, Comparison with normal sensors, advantages, disadvantages, applications.",
        "Module 4": "Actuators: Electric, Hydraulic, Pneumatic. Hydraulic Actuators - Linear and Rotary (Gear, Vane motor). Pneumatic Actuators - Bellows, Flapper-nozzle, Diaphragm actuators. Electric actuators - Solenoids, Stepper motors, DC motors, servo motors. Electro-Pneumatic actuators."
    },

    "ARM Architecture and Programming": {
        "Module 1": "Embedded C: Fixed-width integers in C99, boolean, bit manipulation, memory mapped IO with pointers, structures, bit fields, unions. ARM Cortex-M: Registers, pipelining, memory model, bit banding, instruction format and operands.",
        "Module 2": "ARM assembly programming: Load/store, memory operations, C to ASM conversion (assignments, pointers, subscripts, structures), stack ops, arithmetic, bit manipulation, shifts, bit-field instructions.",
        "Module 3": "Control structures in ASM: instruction sequencing, conditional branches, if-then(-else), loops, function implementation – call/return, register usage, parameter passing, temporary variables, register preservation.",
        "Module 4": "IO programming in assembly: interrupts, exceptions, thread/handler modes, exception handling, latency, priorities, synchronization, buffering, polling, interrupt driven IO, DMA. System initialization: Memory layout, vector table, run-time setup, SysTick."
    },

    "High Speed Digital Design": {
        "Module 1": "Fundamentals: Frequency/time relation, Lumped vs Distributed systems, Reactance types, mutual capacitance/inductance, crosstalk.",
        "Module 2": "High speed logic gates: power dissipation (CMOS), speed, packaging. Measurement techniques: probe effects, shielding, metastability, clock observation.",
        "Module 3": "Transmission lines: point-to-point issues, EMI, crosstalk, lossless/distortionless lines, skin/proximity effects, source/load impedance, terminations, AC biasing, resistor selection.",
        "Module 4": "Power systems: voltage references, bypass capacitors, clock distribution: skew, impedance, duty cycle, oscillators, jitter."
    },

    "Estimation and Detection": {
        "Module 1": "Statistical Estimation I: Estimation basics, MVU estimation, Cramer-Rao Bound, linear models, BLUE, applications.",
        "Module 2": "Statistical Estimation II: MLE, least squares, Bayesian approach, MMSE estimation, applications.",
        "Module 3": "Statistical Detection I: Detection basics, hypothesis testing, Neyman-Pearson, likelihood ratio, ROC, Bayesian detection, Bayes risk, multiple hypothesis.",
        "Module 4": "Statistical Detection II: Detection of deterministic signals (matched filter), random signals, estimator-correlator, linear model applications."
    },

    "ARM Architecture, Programming and Interfacing": {
        "Module 1": "Embedded C: Fixed-width integer types in C99, booleans, bitwise memory/IO access, pointers, structures, packed/bit fields, unions. ARM Cortex-M overview: Registers, pipelining, memory model, bit banding, instruction set.",
        "Module 2": "Assembly programming: Load/store operations, memory addressing, pointer/array/structure translation to ASM, stack, arithmetic/logic/shift instructions.",
        "Module 3": "Assembly control: Instruction sequencing, conditional branching, conditionals and loops, function call/return, register handling.",
        "Module 4": "IO programming: Interrupts, exception handling, latency, priorities, synchronization, buffers, queues, polling, interrupt-driven IO, DMA. System init: Memory map, vector table, SysTick."
    },
    "Advanced Communication Theory" :{
    "Module 1": "Entropy, properties, joint and conditional entropy, mutual information, source coding, Huffman coding, channel capacity (BSC & BEC), Shannon’s channel coding theorem.",
    "Module 2": "AWGN channel capacity, Shannon-Hartley theorem, linear block codes, generator and parity-check matrices, convolutional codes, Viterbi algorithm.",
    "Module 3": "Wireless communication intro, PAN, WiMAX, cellular design fundamentals, handoff, trunking, frequency reuse, FDMA, TDMA, CDMA, OFDMA.",
    "Module 4": "Path loss and shadowing, fading models, OFDM basics, diversity (Alamouti), equalization (linear, nonlinear, MMSE)."
},

"Microwaves and Antennas" : {
    "Module 1": "Microwave spectrum, cavity resonators, hybrid circuits, S-parameters, microwave semiconductor devices (Tunnel, Gunn diode).",
    "Module 2": "Microwave tubes: Klystron, TWT, Magnetron. Measurements of power, frequency, VSWR, Network Analyzer, Anechoic chamber.",
    "Module 3": "Antenna basics, parameters, reciprocity, Helmholtz, directivity and radiation resistance of dipoles.",
    "Module 4": "Antenna arrays, phased arrays, log-periodic, helical, patch, horn, parabolic dish, inverted F antenna."
},

"Computer Networks":{
    "Module 1": "Network basics, switching, queueing models, Little’s theorem, protocol layering, TCP/IP.",
    "Module 2": "Application layer: HTTP, DNS, SMTP. Transport: UDP, TCP, ARQ protocols, congestion control, QoS.",
    "Module 3": "Network layer: Routing, IPv4/IPv6, ARP, CIDR, ICMP, RIP, OSPF, BGP.",
    "Module 4": "Link layer: Error detection (CRC), MAC, Ethernet, CSMA/CD, CSMA/CA, wireless networks, physical media."
},

"Digital Image Processing" :{
    "Module 1": "Image representation, types, pixel relationships, color models, sampling and quantization.",
    "Module 2": "2D transforms: DFT, DCT, Haar, SVD. Compression: JPEG, transform coding.",
    "Module 3": "Image enhancement: Spatial and frequency domain filters, histogram methods.",
    "Module 4": "Restoration: Wiener filtering. Segmentation: Thresholding, edge-based, Hough Transform."
},

"Secure Communication":{
    "Module 1": "Introduction to security, OSI security architecture, symmetric cryptography, Hill & transposition ciphers, finite fields.",
    "Module 2": "Block ciphers: DES, AES, Feistel structure.",
    "Module 3": "Public key cryptography: RSA, Euler’s theorem, key management.",
    "Module 4": "Authentication: MACs, hash functions, digital signatures."
 },

    "NANOELECTRONICS": {
        "Module 1": "Introduction to Nano electronics, Review of MOSFETs, Band diagram operation, threshold voltage, current, MOSFET parameters. Challenges of sub-100 nm MOSFETs, Technological and physical limits, characteristic lengths, scaling and short channel effects, tunneling, power density, dopant concentration, hot electron effects, DIBL, channel length modulation, high-K gate dielectrics and EOT.",
        "Module 2": "Novel MOS devices: SOI (FD, PD), Double Gate MOSFETs, FinFETs, Nanowires. Multi-gate physics, performance optimization with Fins, orientation, Gate stack, gate patterning, Threshold voltage, Work function, Mobility and Strain engineering.",
        "Module 3": "Quantum transport: Atomistic view, quantum of conductance, Coulomb blockade, Schrodinger equation, Band structures, Subbands, Density of states, Ballistic and diffusive transport, Landauer formula.",
        "Module 4": "Applications: Tunneling, MODFET, Resonant tunneling, SET, Coulomb staircase, Spintronics, GMR, TMR, Spin transistor."
    },
    "OPTICAL FIBER COMMUNICATION": {
        "Module 1": "Fiber structure, materials, communication system block diagram, waveguides, types of fibers, attenuation, losses, dispersion types, nonlinear effects.",
        "Module 2": "Fiber fabrication: double crucible, MCVD, OVPO. Fiber cables and splicing, connectors, OTDR, measurements.",
        "Module 3": "Optical sources: LEDs, LDs, detectors (PIN, APD), Receivers, noise analysis, BER, EDFA, Raman amplifiers.",
        "Module 4": "Multiplexing techniques (WDM, OTDM, CDMA etc.), network components, SONET/SDH, LiFi, VLC, sensors."
    },
    "OPTIMIZATION TECHNIQUES": {
        "Module 1": "Optimization introduction, classification, calculus review, stationary points, convexity, LP problems, simplex, duality.",
        "Module 2": "Unconstrained optimization: Fibonacci, Golden Section, Hooke-Jeeves, Newton's method.",
        "Module 3": "Constrained optimization: penalty, barrier methods, Lagrangian, KKT conditions, steepest descent.",
        "Module 4": "Modern optimization: Genetic Algorithms, Simulated Annealing, PSO, Ant Colony. Use of MATLAB/Scilab."
    },
    "IMAGE PROCESSING APPLICATIONS": {
        "Module 1": "Digital image fundamentals, pixel relationships, brightness, contrast, color models, 2D sampling.",
        "Module 2": "2D transforms: DFT, Walsh, Hadamard, Haar, DCT, KL, SVD. Image compression: JPEG, transform coding.",
        "Module 3": "Image enhancement: Spatial and Frequency domain methods, filtering techniques.",
        "Module 4": "Image restoration: inverse filtering, Wiener, constrained filtering. Segmentation: thresholding, clustering, edge, Hough Transform."
    },
    "VLSI CIRCUIT DESIGN": {
        "Module 1": "VLSI Design methodologies, ASIC types, SoCs, FPGA, Top-down/Bottom-up, logical/physical design.",
        "Module 2": "Fabrication: Crystal growth, epitaxy, oxidation, diffusion, lithography, etching, deposition. MOSFET fabrication, CMOS inverter.",
        "Module 3": "CMOS logic design, layout rules, stick diagrams, NAND/NOR layout.",
        "Module 4": "Pass transistor logic, Dynamic logic (Domino, NP Domino), sequential logic, ROM/RAM design."
    },
    "ENTERTAINMENT ELECTRONICS": {
        "Module 1": "Analog TV, PAL/NTSC, digital streaming (MPEG, PSI), Set-top boxes, synchronization.",
        "Module 2": "DVB (S, C, T, H), DVB standards, DAB, Physical layers, SFN, data broadcasting.",
        "Module 3": "HD Video/Audio, video compression (DCT, MPEG-4), blanking intervals, resolution comparison.",
        "Module 4": "Display technologies: CRT, Plasma, LCD, LED, OLED, Field Emission Displays, Holographic, VR, AR."
    },
    "COMPUTER NETWORKS": {
        "Module 1": "Network basics, switching, transmission modes, OSI model, delays, losses, network categories.",
        "Module 2": "TCP/IP: HTTP, SMTP, MIME, POP3, DNS, ARQ protocols, TCP segment, RTT, flow and congestion control.",
        "Module 3": "Network layer: IPv4/IPv6, Routing (static, dynamic), RIP, OSPF, BGP, ARP, RARP, ICMP, CIDR, security.",
        "Module 4": "Link layer: Error detection, CRC, ALOHA, CSMA, MAC, Ethernet, IEEE 802.11, physical media types."
    },
    "BIOMEDICAL ENGINEERING": {
        "Module 1": "Intro, bio-potentials (ECG, EEG, EMG, etc.), electrode theory, types, amplifiers (instrumentation, isolation, chopper).",
        "Module 2": "Heart & cardiovascular system, ECG machine, Nervous system, EEG, EMG, stimulation, myoelectric control.",
        "Module 3": "Clinical lab instruments: oxymeters, counters, spectrophotometers. Therapeutic equipment: pacemakers, defibrillators, dialysis, telemetry.",
        "Module 4": "Medical imaging: X-ray, CT, Ultrasound (A, B, M scans), real-time imaging, MRI, NMR basics."
    },
    "ADVANCED MOBILE COMMUNICATION": {
        "Module 1": "Evolution from 1G to 5G: 1G analog, 2G digital (GSM, CDMA), 2.5G (GPRS), 2.75G (EDGE), 3G (UMTS, W-CDMA, HSPA), 4G (LTE, VoLTE, OFDM, MIMO), LTE Advanced Pro, 5G and enhancements.",
        "Module 2": "5G Basics and applications: eMBB, URLLC, MMTC, D2D, V2X, spectrum sharing, mmWave communication, carrier aggregation, small cells, dual connectivity.",
        "Module 3": "5G Network: New Radio, standalone vs non-standalone, NOMA, massive MIMO, beamforming, SDAP, open RAN, MEC, SDN, NFV, network slicing, restful API, private networks.",
        "Module 4": "Challenges and deployment: 5G in developed and developing countries, backhaul, spectrum, rural connectivity (BharatNet, TVWS, FSO), non-terrestrial solutions (LEO, UAV)."
    },
    "DEEP LEARNING": {
        "Module 1": "Review of ANN, CNN architecture, convolution/pooling layers, feature & weight visualization, t-SNE, AlexNet, VGG, ResNet, GoogLeNet (for practicals).",
        "Module 2": "Loss and activation functions, training CNNs, optimization (SGD, Adam etc.), regularization (dropout, L1/L2), transfer learning, feature extraction.",
        "Module 3": "Sequence models: RNNs, BPTT, vanishing gradients, LSTM, GRU architectures and training.",
        "Module 4": "Generative models: autoencoders, VAEs, GANs, transformers, attention, embeddings, BERT, GPT. (No detailed math needed.)"
    },
    "ROBOTICS AND AUTOMATION": {
        "Module 1": "Basics: robotics vs automation, applications, anatomy, configurations (SCARA, articulated, etc.), workspace, DOF, kinematics, DH notation.",
        "Module 2": "Control systems: open/closed loop, motion types, trajectory planning, controllers (P, PI, PID, etc.).",
        "Module 3": "Actuation & feedback: sensors (encoders, potentiometers), actuators (DC, stepper, servo, hydraulic, pneumatic), grippers, power transmission (gears, chains).",
        "Module 4": "Industrial applications: pick & place, palletizing, die casting, robot cell layouts, safety, control, cycle time analysis."
    },
    "ERROR CONTROL CODING": {
        "Module 1": "Groups, Rings, Finite Fields, Polynomial rings, conjugate/minimal polynomials, vector spaces.",
        "Module 2": "Error control coding: Hamming distance, code rate, LBC, ML decoding, bounds (Singleton, Hamming, etc.), MDS codes.",
        "Module 3": "Cyclic codes: polynomial and matrix view, encoding/decoding, Hamming codes, BCH, Reed-Solomon codes.",
        "Module 4": "Convolution codes, Turbo codes, LDPC (Tanner graphs), message passing, Polar codes in 5G."
    },
    "ADVANCED DIGITAL SIGNAL PROCESSING": {
        "Module 1": "Multirate systems: up/down sampling, polyphase decomposition, filter banks, QMF, perfect reconstruction filters.",
        "Module 2": "Wavelet transforms: STFT, CWT, DWT, Haar wavelet, orthonormal wavelets, time-frequency tradeoffs.",
        "Module 3": "Power spectrum estimation: parametric (Yule-Walker), non-parametric (periodogram), autocorrelation relationships.",
        "Module 4": "Linear prediction filters, adaptive filters (Wiener, LMS), equalization, noise cancellation, LPC of speech."
    },
    "CRYPTOGRAPHY": {
        "Module 1": "Basics: block/stream ciphers, complexity classes (P, NP), number theory (congruences, Euler's, Fermat's), primitive roots.",
        "Module 2": "Algebraic structures, symmetric ciphers (Affine, Hill), DES, AES.",
        "Module 3": "Public key cryptography: RSA, DSA, Knapsack, ECC, zero-knowledge proofs.",
        "Module 4": "Cryptanalysis: primality tests, rho method, linear/differential attacks, factoring (Dixon, Quadratic Sieve)."
    },
    
    "DEEP LEARNING TECHNIQUES": {
        "Module 1": "ANN, Perceptrons, CNNs: convolution, architecture, pooling, FC layers, feature/weight visualization, t-SNE.",
        "Module 2": "Loss functions, activation functions, backpropagation, optimizers (SGD, Adam, etc.), regularization, transfer learning, AlexNet, VGG, ResNet, GoogLeNet (case studies).",
        "Module 3": "RNNs, LSTM, GRU: architectures, training, BPTT, vanishing/exploding gradients.",
        "Module 4": "Generative models: MLE, autoencoders, VAEs, GANs, Transformer models, word embeddings, BERT, GPT (no math)."
    },
    "SATELLITE AND RADAR COMMUNICATION": {
        "Module 1": "Satellite orbits, Kepler's laws, locating satellite, orbital perturbations, launching, communication subsystem.",
        "Module 2": "Satellite link design, uplink power control, interference, rain attenuation, system availability, frequency reuse.",
        "Module 3": "Radar basics: equation, block diagram, performance, CW radar, FM-CW radar, Doppler, MTI vs Pulse Doppler.",
        "Module 4": "Tracking radar techniques."
    },
    "INTERNET OF THINGS": {
        "Module 1": "IoT definitions, architecture, M2M, middleware, service-oriented view, application domains (e.g., smart grid, home, agri, cities).",
        "Module 2": "Addressing (RFID, IPv4/IPv6), sensors, actuators, hardware/software overview, connectivity characteristics.",
        "Module 3": "IoT communication tech: Zigbee, BLE, WiFi, GSM/3G/4G, LoRa, 6LoWPAN, NB-IoT, Sigfox.",
        "Module 4": "IoT data management: storage, flash types, cloud platforms (Azure, Google, IBM), analytics, cloud service models."
    },
    "REAL TIME OPERATING SYSTEM": {
        "Module 1": "Real-time systems overview, RTOS types, task scheduling, timing constraints, IPC (shared memory, message passing), RTOS setup on microcontrollers.",
        "Module 2": "Scheduling: RM, DM, EDF, round robin. Services: memory/I/O management. Middleware. Case studies: automotive, medical.",
        "Module 3": "Inter-task communication: semaphores, mutexes, event flags. Priority inversion and inheritance. Practical synchronization.",
        "Module 4": "RT system design: WCET analysis, modular/time/event-triggered design, fault tolerance, redundancy, system case studies."
    },
    "MIXED SIGNAL CIRCUITS": {
        "Module 1": "CMOS amplifiers: CS, CG, CD, degeneration, cascoded and folded cascode. Small signal models.",
        "Module 2": "MOS current mirrors: simple, cascode. Differential amplifiers: resistive, current source/mirror loads, telescopic cascode.",
        "Module 3": "Op-Amps: two-stage, Miller compensation, bandgap references, supply/temperature independent biasing.",
        "Module 4": "DAC/ADC architectures: resistor string, R-2R, current steering, charge scaling, flash, SAR, pipeline, oversampling."
    },
    "SPEECH AND AUDIO PROCESSING": {
        "Module 1": "Speech production: source-filter model, spectrogram, phonetics, pitch/formants.",
        "Module 2": "Short-time analysis, windowing, STFT, energy/ZCR, filter banks, MFCC, LPC, pitch estimation.",
        "Module 3": "Speech enhancement: spectral subtraction, harmonic filtering. Speaker recognition/identification, ML models, language ID.",
        "Module 4": "Audio perception: anatomy of hearing, filter banks, critical bands, masking, perception models."
    },
    "MICROWAVE DEVICES & CIRCUITS": {
        "Module 1": "Limitations of solid-state devices, IMPATT/TRAPATT diodes, MESFET, amplifier/oscillator design.",
        "Module 2": "Microwave network analysis: S/Y/Z matrices, signal flow graphs, matching (stub, lumped, QWT).",
        "Module 3": "Microwave filters: periodic structures, image/insertion loss method, transformation.",
        "Module 4": "MICs: hybrid vs monolithic, transmission lines (stripline, microstrip), distributed/lumped IC elements."
    },
    "MIXED SIGNAL CIRCUIT DESIGN": {
        "Module 1": "Same as Mixed Signal Circuits Module 1 (CMOS amplifiers, small signal models, cascode, folded).",
        "Module 2": "Same as Mixed Signal Circuits Module 2 (Differential amps, mirrors).",
        "Module 3": "Same as Mixed Signal Circuits Module 3 (Opamps, Miller, Bandgap, DAC/ADC specs).",
        "Module 4": "Same as Mixed Signal Circuits Module 4 (DAC: R-2R, SAR; ADC: Flash, Pipeline, Oversampling)."
    },
    "OPTICAL COMMUNICATION": {
        "Module 1": "Optical fiber basics, structure, advantages, waveguides (NA, V-number, step/graded index), modes, attenuation, dispersion types, nonlinearities, solitons.",
        "Module 2": "Fiber fabrication: double crucible, OVPO, MCVD. Fiber cables, splices, connectors, OTDR, MZ interferometer, attenuation/dispersion measurements.",
        "Module 3": "Optical sources: LEDs, LDs, modulation, coupling. Detectors: PIN, APD, quantum efficiency. Optical receivers, SNR, BER. Amplifiers: EDFA, Raman.",
        "Module 4": "Multiplexing: OTDM, OFDM, WDM, CDMA. SONET/SDH networks. Free space optics, LiFi, VLC, fiber applications in sensors and entertainment."
    },
    "DIGITAL IMAGE PROCESSING": {
        "Module 1": "Image fundamentals: types, pixel relationships, color models (RGB, CMY, HIS), Mach band, brightness, hue, quantization and sampling.",
        "Module 2": "Transforms: DFT, DCT, KL, SVD, Walsh, Hadamard, Haar. Compression: lossy/lossless, transform coding, JPEG standard.",
        "Module 3": "Enhancement: spatial methods (histogram, filters), frequency methods (low/high-pass, homomorphic). Restoration: inverse/Wiener filtering.",
        "Module 4": "Segmentation: thresholding, edge, region, clustering, Hough transform."
    },
    "OPTIMIZATION TECHNIQUES": {
        "Module 1": "Intro to optimization, stationary points, convexity, local/global optima, LP problem formulation, simplex method, duality, LPP applications.",
        "Module 2": "Unconstrained optimization: direct search (Fibonacci, golden section), Hooke-Jeeves, gradient methods (Newton’s).",
        "Module 3": "Constrained optimization: penalty/barrier methods, Lagrangian, KKT conditions, constrained steepest descent.",
        "Module 4": "Modern techniques: genetic algorithms, simulated annealing, PSO, ant colony optimization. MATLAB/Scilab applications."
    }




}










# Add 2024 subjects and modules here

# Combined scheme dictionary

scheme = {
    "2019": syllabus_2019,
    "2024": syllabus_2024
}

import streamlit as st

st.set_page_config(page_title="KTU Study Assistant", layout="centered")

# 💅 Apply colorful custom styles
st.markdown("""
    <style>
        /* 🎨 Main app background (gradient) */
        html, body, [data-testid="stAppViewContainer"] {
            background: linear-gradient(135deg, #f6d365 0%, #fda085 100%);
            color: #000000;
            font-family: 'Comic Sans MS', monospace;
        }

        /* 📚 Sidebar (cool blue with shadow) */
        [data-testid="stSidebar"] {
            background: linear-gradient(180deg, #74ebd5 0%, #ACB6E5 100%);
            color: #000000;
            box-shadow: 2px 0 8px rgba(0,0,0,0.2);
        }

        /* 📘 Title with gradient text */
        h1 {
            font-size: 2.8em;
            font-weight: 800;
            background: linear-gradient(to right, #6a11cb 0%, #2575fc 100%);;
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }

        /* 🌈 Selectboxes and inputs */
        input, textarea, select {
            background-color: #ffffff !important;
            color: #000000 !important;
            border-radius: 12px !important;
            border: 2px solid #ff6f91 !important;
            padding: 6px 10px !important;
        }

        /* 🔴 Buttons */
        button[kind="primary"] {
            background-color: #ff4b5c !important;
            color: white !important;
            font-weight: bold;
            border-radius: 10px !important;
            border: none;
            padding: 8px 16px !important;
            transition: all 0.3s ease-in-out;
        }

        button[kind="primary"]:hover {
            background-color: #ff6f91 !important;
            transform: scale(1.05);
        }

        /* 🟦 Chat input box */
        [data-testid="stChatInput"] textarea {
            background-color: #fff8dc !important;
            border-radius: 10px !important;
            border: 2px dashed #ff9ff3 !important;
        }

        /* 💡 Alert/info box */
        .stAlert {
            background-color: #dfe6e9 !important;
            border-left: 6px solid #00cec9 !important;
        }

        /* 📝 Markdown text */
       
 
    .stMarkdown p {
        font-size: 1.15rem;
        color: #1a1a1a;
        font-weight: 500;
        line-height: 1.6;
        margin: 0;
        padding: 0;
    }
            
    
    </style>
""", unsafe_allow_html=True)









st.title("📘 KTU ECE Study Assistant")
st.markdown("Select a scheme, subject and module, and ask anything from that portion!")

with st.sidebar:
    st.markdown("THIS IS MY DAY")
    selected_year = st.selectbox("Select KTU Scheme", options=["2019", "2024"])
    st.markdown(f"📚 Syllabus: KTU ECE {selected_year} Scheme")
    st.markdown("👩‍🏫 Ask short answers, derivations, MCQs")
    feature = st.selectbox("Choose Feature", [
        "Chatbot Tutor 🤖",
        "Exam Countdown ",
        "SGPA & CGPA Calculator 🎓",
        "PYQ Finder 📂",
        "Quiz Generator 🧠",
        "Notes Summarizer 📝",
        "Study Tracker 📊",
        "Upload & View PDFs 📑",
        "Progress Chart 📈"
    ])

# 🔹 Feature 1: Chatbot
if feature == "Chatbot Tutor 🤖":
    selected_scheme = scheme[selected_year]
    if selected_scheme:
        subject = st.selectbox("Select Subject", list(selected_scheme.keys()), key="chat_subject")
        module = st.selectbox("Select Module", list(selected_scheme[subject].keys()), key="chat_module")

        st.info(f"📘 {subject} - {module}: {selected_scheme[subject][module]}")

        if "chat_history" not in st.session_state:
            st.session_state.chat_history = []

        for role, msg in st.session_state.chat_history:
            with st.chat_message(role):
                st.markdown(msg)

        user_input = st.chat_input(f"Ask from {subject} - {module}...")
        if user_input:
            st.chat_message("user").markdown(user_input)
            st.session_state.chat_history.append(("user", user_input))

            system_prompt = f"""
You are a helpful AI tutor trained in KTU ECE {selected_year} scheme.
The student is asking about **{subject} - {module}**.
Module topics: {selected_scheme[subject][module]}
Answer clearly, concisely, and in KTU exam style with formulas or short examples if needed.
"""
            full_prompt = f"{system_prompt}\n\nUser: {user_input}"
            response = model.generate_content(full_prompt)
            reply = response.text

            st.chat_message("assistant").markdown(reply)
            st.session_state.chat_history.append(("assistant", reply))
    else:
        st.warning("No syllabus added for the selected scheme yet. Please update the 2024 scheme.")

# 🔹 Feature 2: Exam Countdown

elif feature == "Exam Countdown ":
    st.subheader("📅 Exam Countdown & Study Plan")

    today = datetime.date.today()
    exam_date = st.date_input("📅 Select your exam date", value=today)

    if exam_date:
        days_left = (exam_date - today).days
        st.markdown(f"📆 Selected Exam Date: **{exam_date.strftime('%B %d, %Y')}**")

        if days_left > 0:
            st.success(f"⏳ {days_left} days left until your exam!")

            subjects = st.text_area("📚 Enter subjects (comma separated)", "LIC, DSP, VLSI")

            if st.button("Generate Study Plan"):
                subject_list = [s.strip() for s in subjects.split(",") if s.strip()]
                if subject_list:
                    per_subject = max(1, days_left // len(subject_list))
                    st.markdown("### 📘 Suggested Study Plan:")
                    for i, subject in enumerate(subject_list):
                        start_day = i * per_subject + 1
                        end_day = min((i + 1) * per_subject, days_left)
                        st.write(f"• **{subject}**: Day {start_day} to Day {end_day}")
                else:
                    st.warning("⚠️ Please enter at least one subject.")
        elif days_left == 0:
            st.warning("📝 Today is your exam! Best wishes!")
        else:
            st.error("📅 The exam date is in the past. Please choose a future date.")



# 🔹 Feature 3: SGPA & CGPA Calculator
elif feature == "SGPA & CGPA Calculator 🎓":
    st.subheader("🧮 SGPA + CGPA Calculator")

    grade_points = {
        "S": 10.0, "A+": 9.0, "A": 8.5, "B+": 8.0,
        "B": 7.5, "C+": 7.0, "C": 6.5, "D": 6.0,
        "P": 5.5, "F": 0.0, "FE": 0.0, "I": 0.0
    }

    sem = st.selectbox("Select Semester", ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"])
    num_subjects = st.number_input("Number of subjects", min_value=1, max_value=12, step=1)

    subject_data = []
    for i in range(num_subjects):
        col1, col2, col3 = st.columns(3)
        with col1:
            subject = st.text_input(f"Subject {i+1}", key=f"sub_{i}")
        with col2:
            credits = st.number_input("Credits", min_value=1, max_value=5, key=f"cred_{i}")
        with col3:
            grade = st.selectbox("Grade", list(grade_points.keys()), key=f"grade_{i}")
        subject_data.append((subject, credits, grade))

    if st.button("Calculate SGPA"):
        total_credits = sum([c for _, c, _ in subject_data])
        total_points = sum([credits * grade_points[grade] for _, credits, grade in subject_data])
        sgpa = total_points / total_credits if total_credits > 0 else 0
        st.success(f"🎯 SGPA for {sem}: {sgpa:.2f}")

        if "grade_card" not in st.session_state:
            st.session_state.grade_card = {}
        st.session_state.grade_card[sem] = {
            "SGPA": round(sgpa, 2),
            "Subjects": subject_data
        }

    if "grade_card" in st.session_state and st.session_state.grade_card:
        st.markdown("---")
        st.subheader("📁 Saved Semesters")
        total_sgpa = 0
        count = 0
        for sem_name, data in st.session_state.grade_card.items():
            st.markdown(f"### {sem_name} – SGPA: {data['SGPA']}")
            for sub_name, credit, grade in data["Subjects"]:
                st.write(f"• {sub_name} ({credit} credits): Grade {grade}")
            total_sgpa += data["SGPA"]
            count += 1

        if count > 0:
            cgpa = total_sgpa / count
            st.success(f"🎓 CGPA across {count} semesters: {cgpa:.2f}")

            if cgpa >= 8.0:
                classification = "🎓 First Class with Distinction"
            elif cgpa >= 6.5:
                classification = "✅ First Class"
            else:
                classification = "⚠️ Second Class or Below"

            percent = 10 * cgpa - 2.5
            st.info(f"📌 Degree Classification: {classification}")
            st.write(f"📈 Equivalent Percentage: {percent:.2f}%")
elif feature == "PYQ Finder 📂":
    st.subheader("📂 PYQ Search (Coming Soon)")
    st.info("You will be able to search previous year questions by subject or keyword.")

# ✅ Quiz Generator with Gemini
elif feature == "Quiz Generator 🧠":
    st.subheader("🧠 Quiz Generator using Gemini")
    
    subject = st.selectbox("Select Subject", list(syllabus_2019.keys()), key="quiz_subject")
    module = st.selectbox("Select Module", list(syllabus_2019[subject].keys()), key="quiz_module")
    num_questions = st.slider("How many questions?", 1, 10, 5, key="num_qs")
    quiz_type = st.radio("Select Quiz Type", ["Theoretical", "Numerical"], horizontal=True)

    if st.button("Generate Quiz using Gemini"):
        type_text = (
            "conceptual and theory-based multiple choice questions with clear, well-formatted options."
            if quiz_type == "Theoretical"
            else "numerical multiple choice questions with correct calculations. Ensure only one correct answer per question and that the correct option is present."
        )

        prompt = f"""
Generate {num_questions} {type_text} from the topic: {syllabus_2019[subject][module]}.
Format each question like:

Q1. Question text?
A. Option A
B. Option B
C. Option C
D. Option D
Answer: B
        """

        result = model.generate_content(prompt)
        raw_text = result.text
        st.session_state.raw_quiz_text = raw_text

        import re
        question_blocks = re.findall(r"(Q\d+\..*?)(?=\nQ\d+\.|\Z)", raw_text, re.DOTALL)
        parsed_quiz = []

        for block in question_blocks:
            lines = block.strip().split("\n")
            q_line = lines[0].strip()
            q_text = q_line[q_line.find('.')+1:].strip()
            options = {}
            for line in lines[1:]:
                line = line.strip()
                if re.match(r"[A-D]\.", line):
                    label = line[0]
                    content = line[2:].strip()
                    options[label] = content
                elif "Answer:" in line:
                    correct = line.split("Answer:")[-1].strip()
            if len(options) == 4 and correct in options:
                parsed_quiz.append({
                    "question": q_text,
                    "options": options,
                    "answer": correct
                })

        if not parsed_quiz:
            st.error("⚠️ Could not parse quiz. Please try again.")
        else:
            st.session_state.quiz_data = parsed_quiz
            st.session_state.quiz_answers = {}
            st.success("✅ Quiz Generated Successfully!")

    # ✅ Show Quiz Interface
    if "quiz_data" in st.session_state:
        st.markdown("### 📝 Attempt the Quiz")
        for idx, q in enumerate(st.session_state.quiz_data):
            st.markdown(f"**Q{idx+1}. {q['question']}**")
            options = list(q['options'].values())
            user_choice = st.radio(
                f"Select answer for Q{idx+1}", options, key=f"user_ans_{idx}"
            )
            st.session_state.quiz_answers[idx] = user_choice

        if st.button("Submit Quiz"):
            correct = 0
            total = len(st.session_state.quiz_data)
            for idx, q in enumerate(st.session_state.quiz_data):
                user_ans = st.session_state.quiz_answers.get(idx)
                correct_label = q['answer']
                correct_value = q['options'][correct_label]
                if user_ans == correct_value:
                    st.success(f"✅ Q{idx+1} Correct")
                    correct += 1
                else:
                    st.error(f"❌ Q{idx+1} Incorrect — Correct: {correct_value}")
            st.info(f"🎯 Your Score: {correct}/{total}")

# ✅ Notes Summarizer
elif feature == "Notes Summarizer 📝":
    st.subheader("📝 Notes Summarizer with PYQ Extractor")
    uploaded = st.file_uploader("Upload your class notes or PYQ (PDF or .txt)", type=["pdf", "txt"])

    if uploaded:
        text = ""
        if uploaded.name.endswith(".pdf"):
            import fitz  # PyMuPDF
            with fitz.open(stream=uploaded.read(), filetype="pdf") as doc:
                for page in doc:
                    text += page.get_text()
        else:
            text = uploaded.read().decode("utf-8")

        if text.strip() == "":
            st.warning("⚠️ Could not read any text from the file.")
        else:
            model = genai.GenerativeModel("models/gemini-1.5-flash")

            # Prompt for summary and key questions
            prompt = f"""
You are a helpful assistant for summarizing KTU class notes and question papers.

1. Summarize the following text clearly and concisely in KTU exam-ready style.
2. Identify and list:
   - Important 3-mark questions
   - Important 14-mark questions
3. Try to organize questions module-wise if clear from content.
4. Don't repeat content and make it clean to read.

Content:
{text[:4000]}  # truncating for safety
            """

            try:
                with st.spinner("🧠 Generating summary and extracting key questions..."):
                    response = model.generate_content(prompt)
                    summary = response.text

                st.markdown("### ✅ Summary and Important Questions")
                st.markdown(summary)

                # Optional: save summary to download
                st.download_button("💾 Download Summary", summary, file_name="ktu_summary.txt")
            except Exception as e:
                st.error(f"❌ Gemini failed to process: {e}")

# ✅ Study Tracker
elif feature == "Study Tracker 📊":
    import datetime
    st.subheader("📊 Study Tracker")
    st.markdown("Log your daily study hours")

    # Initialize log
    if "study_log" not in st.session_state:
        st.session_state.study_log = []

    # Input hours
    hours = st.slider("Hours studied today", 0, 12, 2)

    if st.button("Add to Log"):
        today = str(datetime.date.today())
        st.session_state.study_log.append((today, hours))
        st.success(f"✅ Logged {hours} hours for {today}")

    # Show current log
    if st.session_state.study_log:
        st.markdown("### 📅 Study Log")
        for date, h in st.session_state.study_log:
            st.write(f"📌 {date} → {h} hour(s)")

# ✅ Upload & View PDFs
elif feature == "Upload & View PDFs 📑":
    st.subheader("📑 Upload & View PDFs")
    file = st.file_uploader("Upload a PDF")
    if file:
        st.markdown("### PDF Preview")
        with fitz.open(stream=file.read(), filetype="pdf") as doc:
            for i, page in enumerate(doc):
                text = page.get_text()
                st.markdown(f"#### Page {i+1}")
                st.text(text)
                if i > 1:
                    break

