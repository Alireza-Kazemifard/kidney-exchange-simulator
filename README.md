# kidney-exchange-simulator-
A Python simulator for the Kidney Exchange problem using the Top Trading Cycles and Chains (TTCC) algorithm, featuring a novel fairness-based policy.
<div align="center">

# **Kidney Exchange TTCC Algorithm Simulator**

**[Persian Version (نسخه فارسی)](./README_FA.md)**

A Python-based simulator for the **Top Trading Cycles and Chains (TTCC)** algorithm for kidney exchange, based on the seminal paper by Roth, Sönmez, and Ünver (2004). This project includes a novel enhancement: a **"Best-Value" (Rule g)** chain selection policy designed to improve fairness.

</div>

---

### **📖 About The Project**

End-Stage Renal Disease is a life-threatening condition where a kidney transplant is often the best treatment. However, there is a severe shortage of donor kidneys, and many patients with a willing living donor cannot receive a transplant due to biological incompatibilities (e.g., blood type or HLA).

**Kidney Paired Donation (KPD)**, or kidney exchange, offers a solution. Incompatible patient-donor pairs can swap donors to find compatible matches. The **TTCC algorithm** is a sophisticated mechanism that facilitates these exchanges by identifying not only direct swaps (**cycles**) but also more complex **chains** that may involve donating to a waitlist.

This project provides a command-line tool to simulate the entire TTCC process, from adding patient-donor pairs to visualizing the final outcomes.

---

### **✨ Key Features**

- **Full TTCC Implementation:** Accurately simulates the core logic of the Top Trading Cycles and Chains algorithm.
- **Realistic Compatibility Checks:** Models both **ABO blood type** compatibility and **HLA crossmatch** (using Panel Reactive Antibody - PRA and unacceptable antigens).
- **Advanced Preference Generation:** Calculates patient preferences for available kidneys based on a utility function that minimizes transplant rejection risk (considering donor age and HLA mismatch).
- **Multiple Chain Selection Rules:** Implements all standard chain selection policies (e.g., max length, priority-based) described in the literature.
- **💡 Innovative "Best-Value" Rule (Rule g):** Our novel contribution, a score-based policy that balances **efficiency** (number of transplants) with **fairness** (prioritizing hard-to-match patients).
- **Interactive CLI:** An easy-to-use menu-driven interface to manage the simulation.
- **State Management:** Save and load the state of the exchange pool to a JSON file.
- **Graph Visualization:** Generates graphs of the exchange pool at different stages using NetworkX and Matplotlib.

---

### **🚀 Getting Started**

Follow these steps to get the simulator running on your local machine.

#### **Prerequisites**

- Python 3.6 or higher
- Pip package manager

#### **Installation**

1.  **Clone the repository:**
    ```sh
    git clone https://github.com/YOUR_USERNAME/YOUR_REPOSITORY_NAME.git
    cd YOUR_REPOSITORY_NAME
    ```

2.  **Install the required packages:**
    The project requires `networkx` for graph operations and `matplotlib` for visualization.
    ```sh
    pip install networkx matplotlib
    ```

#### **How to Run**

Execute the main script from your terminal:
```sh
python kidney.py
```
You will be greeted with the main menu, from which you can add pairs, generate preferences, run the algorithm, and more.

---

### **⚙️ How It Works**

The simulation is built around two core concepts from the TTCC algorithm:

- **Cycles:** A closed loop of exchanges where a group of patients directly swap kidneys among themselves. For example, `P1 -> K2 -> P2 -> K1 -> P1`. Cycles are self-sufficient and are the highest priority for execution.
- **w-Chains:** An open path of exchanges that begins with a patient-donor pair and ends with a donation to the general waitlist (`w`). For example, `P1 -> K2 -> P2 -> w`. Chains are considered only when no cycles can be formed.

#### **Chain Selection Rules**

When multiple w-chains are available, a policy is needed to choose one. This simulator implements several rules:

- **`a` (Min Length, Remove):** Chooses the shortest chain and finalizes it.
- **`b` (Max Length, Remove):** Chooses the longest chain and finalizes it (greedy approach).
- **`c` (Max Length, Keep):** Chooses the longest chain but keeps its participants active, allowing the chain to potentially grow or become part of a cycle later.
- **`d` (Priority, Remove):** Chooses the first chain containing the highest-priority patient.
- **`e` (Priority, Keep):** Similar to `d`, but keeps participants active.
- **`f` (Hybrid O-type):** A policy that prioritizes removing chains that donate a valuable Type O kidney to the waitlist.

#### **⭐ Rule g: Best-Value (Score-Based) - Our Contribution**

This project introduces an innovative rule designed to provide a more holistic and fair matching outcome. Instead of focusing on a single metric, it assigns a score to each potential chain based on a combination of factors:

- **+10 points** for each patient in the chain (prioritizes efficiency).
- **+5 points** for each Type O patient in the chain (addresses the disadvantage of this blood type).
- **+10 points** for each highly-sensitized (High PRA) patient (prioritizes hard-to-match cases).

The chain with the highest total score is selected. This allows the algorithm to make more nuanced decisions that balance the overall number of transplants with the urgent needs of specific patient groups.

---

### **📂 Project Structure**

```
.
├── kidney.py                # Main Python script containing all the logic.
├── sample.json              # A sample data file with 60 patient-donor pairs.
├── graph_initial_state.png  # Example output graph.
├── README.md                # This file (English).
└── README_FA.md             # This file (Persian).
```

---

### **🙏 Acknowledgments**

This work is a direct implementation and extension of the concepts presented in the following foundational paper:
- Roth, A. E., Sönmez, T., & Ünver, M. U. (2004). Kidney exchange. *The Quarterly Journal of Economics*, 119(2), 457-488.

---

### **📜 License**

This project is licensed under the MIT License. See the `LICENSE` file for details.
