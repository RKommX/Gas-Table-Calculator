## Gas Table Calculator

The calculator uses the following packages:
- math: For accessing some basic calculations and constants
- numpy: Used in oblique_shock.py to interpolate wave angle with turn angle
- scipy: To access scipy.optimize.newton (Newton-Raphson method) which was used in many inverse functions
- streamlit: For building the GUI

To run the GUI,
- Make sure the above packages are installed
- Download and extract the repository
- Navigate to this directory in the terminal, and run `streamlit run GUI.py`