User Guide: Composite Laminate and Beam Analysis Tool
This tool is a MATLAB-based suite designed for the structural analysis of composite laminates and open-section beams, supporting engineers in verifying designs directly from load cases. It implements Classical Lamination Theory (CLT) to compute full ABD matrices and performs Tsai–Wu ply-by-ply failure checks.
The main.m file is used to analyze laminates made of composite layup, whereas the Beam_analysis_unsymmetric_two.m code handles open-section beams. The rest are helper functions used in these two main functions. 
________________________________________
1. Coordinate Systems
To ensure accurate results, inputs must follow these two coordinate conventions:
A. Global/User Coordinate System (Y, Z)
This system is used to define the cross-sectional geometry of the beam.
•	Origin (0,0): Located at the bottom-left of the section.
•	Z-Axis (Horizontal): Positive is to the LEFT.
•	Y-Axis (Vertical): Positive is UP.
•	Web Offset: The z_web_center input allows the web to be placed anywhere along the width, enabling unsymmetric analysis.
B. Internal/Local Coordinate System
Mechanical calculations (bending and shear) are automatically performed relative to the Modulus-Weighted Centroid (Ybar, Zbar).
•	The code converts your user-defined coordinates into distances from the centroid (note the neutral axis may be at any angle due to bending moments acting in different directions) axis to calculate stresses like σxx.
________________________________________
2. Input Requirements: What to Enter
The script follows a modular flow to gather necessary design data:
Phase 1: Geometry & Materials
•	Dimensions: Enter flange widths and the total web height in inches.
•	Material Cards: The tool uses built-in or user-extended material libraries for unidirectional composite properties (Exx, Eyy, Gxy, etc.) and strength allowables.
Phase 2: Laminate Layups
•	Sequence: Enter plies as angle-count pairs (e.g., 45 2 -45 2 0 4) followed by 000 to terminate the input.
•	Thickness: Provide the cured ply thickness (tply), which determines the total laminate thickness used for the A, B, D matrices.
Phase 3: Applied Loads
•	Bending Moments: Myy (about vertical) and Mzz (about horizontal) in in-lb.
•	Shear Forces: Vy (vertical) and Vz (horizontal) in lb.
•	Axial Load: Px in lb (Tension is positive).
________________________________________
3. What Results to Expect
The tool provides both numerical verification and visual insights into the beam's behavior:
A. Modulus-Weighted Properties
•	Calculates the effective stiffness-weighted area (A*) and inertia tensor (Iyy, Izz, Iyz).
•	Determines the Shear Center, which is critical for understanding torsion-bending coupling in open sections.
B. Stress & Force Distributions (Heatmaps)
•	Normal Force Resultant (Nxx): A color-mapped plot showing axial tension and compression across the section.
•	Shear Flow (Nxy): Visualizes the load paths and shear distribution, identifying where the section is most heavily taxed15.
C. Progressive Failure Summary
•	Ply-by-Ply Analysis: Evaluates failure indices at specific user-defined points.
•	Failure Order: Provides a table showing which plies fail first and the corresponding Load Factor (R).
•	R-Values: Strength ratios are plotted for every ply, where values below 1.0 indicate failure.


