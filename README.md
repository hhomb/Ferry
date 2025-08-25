# Design of a Decision Support System for Energy-Efficient Ferry Operation on Lake Constance based on Optimal Control

📄 **Submitted to the European Control Conference (ECC) 2026**  

---

## Overview  
This repository accompanies the paper:  
**"Design of a Decision Support System for Energy-Efficient Ferry Operation on Lake Constance based on Optimal Control"**  

The work addresses the **energy-efficient operation of inland ferries** by combining:  
- 🚢 A dynamic model of the ferry *MS Insel Mainau*  
- 🎯 A shrinking-horizon optimal control framework  
- 🌊 Environmental disturbances (water currents, wind)  
- ⚡ Real-time decision support for ferry crews  

---

## Abstract  
The maritime sector is undergoing a disruptive technological change driven by autonomy, decarbonization, and digitalization. Addressing these factors necessitates a reassessment of inland vessel operations.  

This paper presents the **design and technical development of a decision support system for ferry operation** grounded in a shrinking-horizon optimal control framework. The formulation integrates a **mathematical model of the ferry’s dynamics** and **environmental disturbances**, specifically water currents and wind, which can significantly affect operations.  

Illustrative experiments demonstrate the system’s potential to:  
- Provide **real-time, actionable guidance** (nudges) to the crew  
- Enhance **operational efficiency**  
- Maintain predefined maneuver durations  
- Ensure **computational tractability** for real-time applications  

The findings suggest that optimal control can advance and potentially transform ferry operations on inland waters.  

---

## Repository Structure  
🔒 **Note:** The code and ferry model will be released after acceptance of the paper.  

Planned contents include:  
- `model/` – Dynamic model of the ferry *MS Insel Mainau*  
- `control/` – Optimal control formulations and solver setup  
- `experiments/` – Simulation and evaluation scripts  
- `docs/` – Additional documentation and results  

---

## Video Demonstration  
🎥 A video of the real-world ferry *MS Insel Mainau* operating on Lake Constance is provided here:  
https://youtu.be/i1MjCdbEQyE
