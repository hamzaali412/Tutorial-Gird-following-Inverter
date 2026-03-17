# Tutorial-Gird-following-Inverter
A comprehensive Python tutorial for simulating grid-following inverters. This project implements Differential-Algebraic Equations (DAEs) solved via the Implicit Trapezoidal Rule and Newton Raphson method.  There are three different operational support modes: Constant PQ, Frequency support and Volt-Var 

## 📐 System Configuration

### Network Topology
The system consists of a two-bus network:
* **Slack Bus (Grid):** An infinite bus with a fixed voltage magnitude ($V_{slack}$) and angle ($\theta_{slack}$).
* **Inverter Bus:** The point of common coupling (PCC) where the inverter is connected.

### Physical Interface
The inverter is interfaced with the grid through a **series RL filter**. The output voltage of the inverter, $v_g^{abc}$, is dynamically controlled to manage the power flow across the filter's coupling impedance to the second node.

### Control Architecture
The control system is implemented in the $dq$-reference frame and consists of:
1. **Synchronization:** A Phase-Locked Loop (PLL) to track the grid frequency and angle.
2. **Outer Power Loop:** Manages Active ($P$) and Reactive ($Q$) power setpoints.
3. **Inner Current Loop:** Fast-acting control of $i_d$ and $i_q$.



<div align="center">
  <a href="diagram/main_diagram.png">
    <img 
      src="diagram/main_diagram.png" 
      alt="System Control Diagram" 
      style="background-color: #e6ebf5; 
             padding: 15px; 
             border-radius: 10px; 
             width: 75%; 
             border: 2px solid #5c7cfa;
             box-shadow: 0 0 20px rgba(92, 124, 250, 0.3);"
    >
  </a>
  <p align="center">
    <br>
    <code style="color: #5c7cfa;">Figure 1: Control block diagram and physical system layout.</code>
  </p>
</div>