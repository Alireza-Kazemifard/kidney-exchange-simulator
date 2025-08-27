# Alireza Kazemifard
import math
import json
import os
import textwrap

# ############################################################################
# 1. HLA Universe and Data Structures
# ############################################################################
HLA_UNIVERSE = {
    'A': ['A1', 'A2', 'A3', 'A11', 'A23', 'A24', 'A26','A29', 'A30', 'A31', 'A32', 'A68'],
    'B': ['B7', 'B8','B13','B15', 'B27', 'B35', 'B40', 'B44', 'B51', 'B57', 'B60', 'B62'],
    'DR': ['DR1', 'DR3', 'DR4', 'DR7', 'DR11', 'DR13', 'DR15', 'DR17', 'DR51', 'DR52']
}

class Patient:
    """Represents a patient with a willing donor."""
    def __init__(self, patient_id, blood_type, age, hla_profile, donor_kidney_id, pra=0, unacceptable_antigens=None , wants_waitlist=False):
        self.id = patient_id
        self.blood_type = blood_type
        self.age = age
        self.hla_profile = hla_profile
        self.donor_kidney_id = donor_kidney_id
        self.pra = pra
        self.unacceptable_antigens = unacceptable_antigens if unacceptable_antigens is not None else set()
        self.preferences = []
        self.is_active = True
        self.assignment = None
        self.wants_waitlist = wants_waitlist

    def __repr__(self):
        return f"Patient(id={self.id}, type='{self.blood_type}', age={self.age}, pra={self.pra}%)"

class Kidney:
    """Represents a donor's kidney."""
    def __init__(self, kidney_id, blood_type, age, hla_profile, donor_patient_id):
        self.id = kidney_id
        self.blood_type = blood_type
        self.age = age
        self.hla_profile = hla_profile
        self.donor_patient_id = donor_patient_id

    def __repr__(self):
        return f"Kidney(id={self.id}, type='{self.blood_type}', donor_age={self.age})"

# ############################################################################
# 2. Main Kidney Exchange Class
# ############################################################################
class KidneyExchange:
    """Orchestrates the kidney exchange using the TTCC algorithm."""
    def __init__(self):
        self.patients = {}
        self.kidneys = {}
        self.next_id = 1
        self.verbose = True

    def _log(self, message):
        if self.verbose:
            print(message)

    # ########################################################################
    # Section 2.1: Data Management
    # ########################################################################
    def reset(self):
        self.patients.clear()
        self.kidneys.clear()
        self.next_id = 1
        self._log("\nSystem has been reset.")

    def add_pair(self, patient_data, donor_data, wants_waitlist=False):
        patient_id = f"p{self.next_id}"
        kidney_id = f"k{self.next_id}"
        
        patient = Patient(
            patient_id, patient_data['blood_type'], patient_data['age'], patient_data['hla'], kidney_id,
            pra=patient_data.get('pra', 0),
            unacceptable_antigens=patient_data.get('unacceptable_antigens', set()),
            wants_waitlist=wants_waitlist
        )
        self.patients[patient_id] = patient

        kidney = Kidney(kidney_id, donor_data['blood_type'], donor_data['age'], donor_data['hla'], patient_id)
        self.kidneys[kidney_id] = kidney
        self.next_id += 1
        return patient_id

    def save_state(self, filepath):
        if not self.patients:
            self._log("Nothing to save. The system is empty."); return
        
        patients_data = [{
            'patient_id': p.id, 'blood_type': p.blood_type, 'age': p.age,
            'hla_profile': p.hla_profile, 'donor_kidney_id': p.donor_kidney_id,
            'pra': p.pra, 'unacceptable_antigens': list(p.unacceptable_antigens),
            'wants_waitlist': p.wants_waitlist,
            'assignment' : p.assignment
        } for p in self.patients.values()]

        kidneys_data = [{'kidney_id': k.id, 'blood_type': k.blood_type, 'age': k.age,
                         'hla_profile': k.hla_profile, 'donor_patient_id': k.donor_patient_id}
                        for k in self.kidneys.values()]

        with open(filepath, 'w') as f:
            json.dump({'patients': patients_data, 'kidneys': kidneys_data, 'next_id': self.next_id}, f, indent=4)
        self._log(f"System state saved to {filepath}")

    def load_state(self, filepath):
        if not os.path.exists(filepath):
            self._log(f"Error: File not found at {filepath}"); return
        
        self.reset()
        with open(filepath, 'r') as f:
            data = json.load(f)

        for p_data in data.get('patients', []):
            unacceptable = set(p_data.pop('unacceptable_antigens', []))
            pra_val = p_data.pop('pra', 0)
            p_data.pop('assignment',[])
            p_id = p_data['patient_id']
            self.patients[p_id] = Patient(**p_data, pra=pra_val, unacceptable_antigens=unacceptable)

        for k_data in data.get('kidneys', []):
            self.kidneys[k_data['kidney_id']] = Kidney(**k_data)
            
        self.next_id = data.get('next_id', self.next_id)
        self._log(f"System state loaded from {filepath}")
        self.list_all_pairs()

    # ########################################################################
    # Section 2.2: Preference Generation with Virtual Crossmatch
    # ########################################################################
    def _is_blood_compatible(self, patient_bt, kidney_bt):
        if kidney_bt == 'O': return True
        if kidney_bt == 'A' and patient_bt in ['A', 'AB']: return True
        if kidney_bt == 'B' and patient_bt in ['B', 'AB']: return True
        return kidney_bt == 'AB' and patient_bt == 'AB'

    def _is_crossmatch_negative(self, patient, kidney):
        donor_hlas = set()
        for locus in kidney.hla_profile.values():
            for hla in locus:
                donor_hlas.add(hla)
        return not donor_hlas.intersection(patient.unacceptable_antigens)

    def _calculate_hla_mismatch(self, patient_hla, donor_hla):
        mismatches = 0
        for locus in ['A', 'B', 'DR']:
            patient_antigens = set(patient_hla.get(locus, []))
            for donor_antigen in donor_hla.get(locus, []):
                if donor_antigen not in patient_antigens:
                    mismatches += 1
        return mismatches

    def _calculate_utility(self, patient, kidney):
        a, b = (1.06, 1.12) if patient.age < 60 else (1.05, 1.10)
        hla_mismatch = self._calculate_hla_mismatch(patient.hla_profile, kidney.hla_profile)
        return -math.log(a) * hla_mismatch - (math.log(b) * kidney.age / 10.0)

    def generate_all_preferences(self):
        if not self.patients: self._log("No patients."); return
        self._log("\nGenerating preferences (Blood Type + Virtual Crossmatch)...")
        for patient in self.patients.values():
            utility_scores = []
            for k_id, kidney in self.kidneys.items():
                
                if self._is_blood_compatible(patient.blood_type, kidney.blood_type) and \
                   self._is_crossmatch_negative(patient, kidney):
                    utility = self._calculate_utility(patient, kidney)
                    utility_scores.append((utility, k_id))
            
            utility_scores.sort(key=lambda item: item[0], reverse=True)
            patient.preferences = [kid for score, kid in utility_scores]
            
            if patient.wants_waitlist: patient.preferences.append('w')
            
        self._log("All preferences generated."); self.list_all_pairs()
    
    # ########################################################################
    # Section 2.3: Core TTCC Algorithm Logic
    # ########################################################################
    def _build_pointers(self, active_patients, available_kidneys):
        pointers = {}
        for p_id in active_patients:
            patient = self.patients[p_id]
            for pref in patient.preferences:
                if pref == 'w' or pref in available_kidneys:
                    pointers[p_id] = pref; break
        for k_id in available_kidneys:
            pointers[k_id] = self.kidneys[k_id].donor_patient_id
        return pointers


    def _execute_cycle(self, cycle):
        self._log(f"  Executing Cycle: {' -> '.join(cycle)} -> {cycle[0]}")
        for i, node_id in enumerate(cycle):
            if node_id.startswith('p'):
                patient = self.patients[node_id]
                patient.assignment = cycle[(i + 1) % len(cycle)]
                patient.is_active = False

    def _process_chain(self, chain, rule):
        self._log(f"  Processing Chain: {' -> '.join(chain)}")
        for i, node_id in enumerate(chain[:-1]):
            if node_id.startswith('p'):
                self.patients[node_id].assignment = chain[i + 1]
        
        action = 'remove'
        if rule in ['c', 'e','g']: action = 'keep'
        elif rule == 'f':
            tail_donor_bt = self.kidneys[self.patients[chain[0]].donor_kidney_id].blood_type
            action = 'remove' if tail_donor_bt == 'O' else 'keep'
            self._log(f"    (Hybrid rule 'f': Tail donor is {tail_donor_bt}, chain will be {'removed' if action == 'remove' else 'kept'}.)")
        
        if action == 'remove':
            self._log("    (Chain finalized. Participants are now inactive.)")
            for node_id in chain:
                if node_id.startswith('p'): self.patients[node_id].is_active = False
        else:
            self._log("    (Chain kept. Participants remain active for potential future cycles.)")

    def _filter_and_select_chain(self, chains, rule, priority_list, max_len):
        filtered = [c for c in chains if len([n for n in c if n.startswith('p')]) <= max_len]
        if not filtered: return None

        self._log(f"  Found {len(filtered)} w-chain(s) within length limit.")
        
        if rule == 'g':
            self._log("  Rule 'g' (Best-Value): Scoring chains...")
            def score_chain(chain):
                score = len([n for n in chain if n.startswith('p')]) * 10
                for node in chain:
                    if node.startswith('p'):
                        p = self.patients[node]
                        if p.blood_type == 'O': score += 5
                        if p.pra >= 80: score += 10
                return score
            
            scored_chains = sorted([(score_chain(c), c) for c in filtered], key=lambda x: x[0], reverse=True)
            self._log(f"    Top chain scored {scored_chains[0][0]}: {' -> '.join(scored_chains[0][1])}")
            return scored_chains[0][1]

        filtered.sort(key=lambda x: (-len(x), x[0]))

        if rule == 'a': return min(filtered, key=len)
        if rule in ['b', 'c']: return filtered[0]
        if rule in ['d', 'e', 'f']:
            for p_id in priority_list:
                for chain in filtered:
                    if p_id in chain: return chain
            return None
        return filtered[0]

    def _find_cycles_and_w_chains(self, pointers, active_patients, processed_in_keep=None):
        if processed_in_keep is None:
            processed_in_keep = set()
            
        cycles, chains, visited = [], [], set()
        for p_id in active_patients:
            if p_id in visited:
                continue
                
            path, curr = [], p_id
            while curr in pointers and curr not in path and curr not in processed_in_keep:
                path.append(curr)
                curr = pointers.get(curr)
            
            if curr == 'w':
                chains.append(path + ['w'])
            elif curr in processed_in_keep:
                chains.append(path + [curr])
            elif curr in path:
                cycle = path[path.index(curr):]
                cycles.append(cycle)
   
            visited.update(path)
        
        unique_cycles, seen_nodes = [], set()
        for cycle in cycles:
            has_overlap = False 
            for node in cycle:
                if node in seen_nodes:
                    has_overlap = True
                    break 
            if not has_overlap:
                unique_cycles.append(cycle)
                seen_nodes.update(cycle)
                
        return unique_cycles, chains

    def _unravel_kept_chain(self, start_patient_id):
        path = [start_patient_id]
        curr_p = self.patients[start_patient_id]
        
        while curr_p.assignment and curr_p.assignment != 'w':
            assigned_k_id = curr_p.assignment
            next_p_id = self.kidneys[assigned_k_id].donor_patient_id
            
            path.append(assigned_k_id)
            path.append(next_p_id)
            
            curr_p = self.patients[next_p_id]
        
        if curr_p.assignment == 'w':
            path.append('w')
            
        return path

    def run_ttcc(self, chain_rule='c', max_cycle_len=999, max_chain_len=999):
        for p in self.patients.values():
            p.is_active, p.assignment = True, None
        priority_list = sorted(self.patients.keys(), key=lambda x: int(x[1:]))
        
        initial_pointers = self._build_pointers(set(self.patients.keys()), set(self.kidneys.keys()))
        self._generate_graph('graph_initial_state', initial_pointers, set(self.patients.keys()))
        
        round_num, processed_in_keep, post_cycle_graph_generated, no_cycles_found = 1, set(), False , 0
        while True:
            self._log(f"\n--- Round {round_num} ---")
            
            active_for_cycles = {pid for pid, p in self.patients.items() if p.is_active}
            if not active_for_cycles:
                self._log("No active patients remain. Algorithm finished.")
                break

            assigned_k = {p.assignment for p in self.patients.values() if p.assignment and p.assignment != 'w'}
            available_k = {k_id for k_id, k in self.kidneys.items() if self.patients[k.donor_patient_id].is_active and k_id not in assigned_k}
            
            pointers_for_cycles = self._build_pointers(active_for_cycles, available_k)
            cycles, _ = self._find_cycles_and_w_chains(pointers_for_cycles, active_for_cycles)
            
            filtered_cycles = [c for c in cycles if (len(c) / 2) <= max_cycle_len]
            if filtered_cycles:
                self._log(f"Found {len(filtered_cycles)} cycle(s) to execute.")
                for cycle in filtered_cycles:
                    self._execute_cycle(cycle)
                round_num += 1
                post_cycle_graph_generated = False
                continue
            else :
                no_cycles_found += 1
                
            if not post_cycle_graph_generated and no_cycles_found == 1:
                self._log("--- Drawing Intermediate State Graph (Post-Cycles) ---")
                pointers_post_cycles = self._build_pointers(active_for_cycles, available_k)
                self._generate_graph('graph_post_cycles_state', pointers_post_cycles, active_for_cycles)
                post_cycle_graph_generated = True
            
            self._log("No cycles found. Looking for a chain to execute...")
            
            active_for_chains = {pid for pid in active_for_cycles if pid not in processed_in_keep}
            pointers_for_chains = self._build_pointers(active_for_chains, available_k)
            _, chains = self._find_cycles_and_w_chains(pointers_for_chains, active_for_chains, processed_in_keep)

            expanded_chains = []
            for chain in chains:
                last_node = chain[-1]
                if last_node in processed_in_keep:
                    unraveled_part = self._unravel_kept_chain(last_node)
                    expanded_chains.append(chain[:-1] + unraveled_part)
                else:
                    expanded_chains.append(chain)

            if expanded_chains:
                selected = self._filter_and_select_chain(expanded_chains, chain_rule, priority_list, max_chain_len)
                if not selected:
                    self._log("No selectable w-chains found."); break
                
                self._process_chain(selected, chain_rule)
                
                action = 'remove'
                if chain_rule in ['c', 'e' , 'g']:
                    action = 'keep'
                elif chain_rule == 'f':
                    initiating_patient_id = selected[0]
                    tail_donor_kidney_id = self.patients[initiating_patient_id].donor_kidney_id
                    tail_donor_bt = self.kidneys[tail_donor_kidney_id].blood_type
                    action = 'remove' if tail_donor_bt == 'O' else 'keep'

                patients_in_chain_path = {node for node in selected if node.startswith('p')}
                
                
                remaining_active_after_chain = active_for_cycles - patients_in_chain_path
                if not remaining_active_after_chain:
                    self._log("    (This is the final transaction, overriding action to 'remove'.)")
                    action = 'remove'
                

                if action == 'keep':
                    processed_in_keep.update(patients_in_chain_path)
                else:
                    for p_id in patients_in_chain_path:
                        self.patients[p_id].is_active = False
            else:
                self._log("No cycles or chains found. Algorithm finished.")
                break
            round_num += 1
            
        self._log("\n--- Drawing Final State Graph ---")
        final_unmatched_patients = {pid for pid, p in self.patients.items() if not p.assignment}
        final_available_kidneys = {k_id for k_id, k in self.kidneys.items() if k.donor_patient_id in final_unmatched_patients}
        final_pointers = self._build_pointers(final_unmatched_patients, final_available_kidneys)
        self._generate_graph('graph_final_state', final_pointers, final_unmatched_patients)
        
        for p in self.patients.values():
            if not p.assignment:
                p.assignment = p.donor_kidney_id
        self._log("\n--- Final Results ---")
        self.display_final_results()

    # ########################################################################
    # Section 2.4: Display and Visualization
    # ########################################################################
    def list_all_pairs(self):
        if not self.patients: print("No patient-donor pairs."); return
        print("\nCurrent Patient-Donor Pairs:"); print("-" * 95)
        for p_id in sorted(self.patients.keys(), key=lambda x: int(x[1:])):
            p = self.patients[p_id]
            k = self.kidneys[p.donor_kidney_id]
            
            print(f"Pair {p_id[1:]:<2}: Patient {p.id} (Type {p.blood_type}, Age {p.age}, PRA {p.pra}%) | Donor {k.id} (Type {k.blood_type}, Age {k.age})")
            print(f"  - Unacceptable Antigens: {sorted(list(p.unacceptable_antigens)) if p.unacceptable_antigens else 'None'}")
            
            if p.preferences:
                prefs_str = ", ".join(p.preferences)
                wrapped_prefs = '\n'.join(textwrap.wrap(prefs_str, 86, initial_indent='  - Preferences: [', subsequent_indent=' ' * 17, break_long_words=False))
                print(f"{wrapped_prefs}]")
            
            print("-" * 95)

    def display_final_results(self):
        if not self.patients: self._log("System is empty."); return
        headers = ["Pair", "Patient", "Pt. Age", "Pt. Type", "Donor", "Dn. Age", "Dn. Type", "Outcome"]
        print(f"\n{headers[0]:<5} {headers[1]:<8} {headers[2]:<8} {headers[3]:<8} {headers[4]:<8} {headers[5]:<8} {headers[6]:<8} {headers[7]:<15}")
        print("-" * 80)
        
        transplants, waitlisted, no_exchange = 0, 0, 0
        for p_id in sorted(self.patients.keys(), key=lambda x: int(x[1:])):
            p = self.patients[p_id]
            k = self.kidneys[p.donor_kidney_id]
            outcome = p.assignment
            display_outcome = ""

            if outcome == 'w':
                waitlisted += 1; display_outcome = "Waitlisted"
            elif outcome == p.donor_kidney_id:
                no_exchange += 1; display_outcome = "No Exchange"
            else:
                transplants += 1; display_outcome = outcome
            
            print(f"{p_id[1:]:<5} {p.id:<8} {p.age:<8} {p.blood_type:<8} {k.id:<8} {k.age:<8} {k.blood_type:<8} {str(display_outcome):<15}")
        
        print("-" * 80); self._log(f"Summary: {transplants} successful transplant(s), {waitlisted} patient(s) to waitlist, {no_exchange} pair(s) with no exchange.")

    def _generate_graph(self, filename, pointers, active_nodes):
        try:
            import networkx as nx; import matplotlib.pyplot as plt
        except ImportError: self._log("\n[NetworkX/Matplotlib not found for graph generation.]"); return

        G = nx.DiGraph()
        all_patients = sorted(list(self.patients.keys()), key=lambda x: int(x[1:]))
        all_kidneys = sorted(list(self.kidneys.keys()), key=lambda x: int(x[1:]))
        G.add_nodes_from(all_patients); G.add_nodes_from(all_kidneys); G.add_node('w')

        pos, node_colors, labels = {}, [], {}
        
        p_angles = {p_id: i * 2 * math.pi / len(all_patients) for i, p_id in enumerate(all_patients)} if all_patients else {}
        for p_id, angle in p_angles.items(): pos[p_id] = (math.cos(angle), math.sin(angle))
        k_angles = {k_id: i * 2 * math.pi / len(all_kidneys) for i, k_id in enumerate(all_kidneys)} if all_kidneys else {}
        for k_id, angle in k_angles.items(): pos[k_id] = (0.6 * math.cos(angle), 0.6 * math.sin(angle))
        pos['w'] = (0, 0)
        
        for node in G.nodes():
            if node.startswith('p'):
                node_colors.append('lightgreen' if node in active_nodes else '#cccccc')
                labels[node] = f"{node}\n({self.patients[node].blood_type}, {self.patients[node].pra}%)"
            elif node.startswith('k'):
                node_colors.append('skyblue'); labels[node] = f"{node}\n({self.kidneys[node].blood_type})"
            else:
                node_colors.append('gold'); labels[node] = 'Waitlist'
        
        edges_to_draw = [(u, v) for u, v in pointers.items() if (u.startswith('p') and u in active_nodes) or (u.startswith('k') and self.kidneys[u].donor_patient_id in active_nodes)]
    
        G.add_edges_from(edges_to_draw)

        plt.figure(figsize=(16, 16))
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=2500, edgecolors='black')
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=9)
        nx.draw_networkx_edges(G, pos, edgelist=G.edges(), edge_color='gray', arrows=True, arrowsize=20, node_size=2500)
        plt.title(f"Kidney Exchange Graph - {filename.replace('_', ' ').title()}", fontsize=16)
        plt.tight_layout(); plt.savefig(f"{filename}.png", format="PNG"); plt.close()
        self._log(f"\nGraph saved to '{filename}.png' using Matplotlib.")
##########################################################################
# paper example
##########################################################################

    def _load_paper_example(self):
        """Loads the 12-pair example from the Roth, Sonmez, Unver (2004) paper."""
        self.reset()
        self._log("\n--- Loading 12-Pair Example from TTCC Paper ---")
        p_types = ['AB', 'O', 'A', 'B', 'A', 'O', 'B', 'A', 'O', 'AB', 'O', 'B']
        d_types = ['B', 'A', 'O', 'A', 'O', 'O', 'A', 'O', 'B', 'B', 'O', 'A']
  
        for i in range(1, 13):
            p_id = f"p{i}"
            k_id = f"k{i}"
            patient = Patient(p_id, p_types[i-1], 40, {}, k_id)
            self.patients[p_id] = patient
            self.kidneys[k_id] = Kidney(k_id, d_types[i-1], 40, {}, p_id)
    
        self.next_id = 13
        self.patients['p1'].preferences = ['k9', 'k10', 'k1']
        self.patients['p2'].preferences = ['k11', 'k3', 'k5', 'k6', 'k2']
        self.patients['p3'].preferences = ['k2', 'k4', 'k5', 'k6', 'k7', 'k8', 'k11','k12', 'w']
        self.patients['p4'].preferences = ['k5', 'k9', 'k1', 'k8', 'k10','k3' ,'k6','w']
        self.patients['p5'].preferences = ['k3', 'k7', 'k11', 'k4', 'k5']
        self.patients['p6'].preferences = ['k3', 'k5', 'k8', 'k6']
        self.patients['p7'].preferences = ['k6', 'k11', 'k1', 'k3','k9','k10', 'k1', 'w']
        self.patients['p8'].preferences = ['k6', 'k4', 'k11', 'k2', 'k3', 'k8']
        self.patients['p9'].preferences = ['k3', 'k11', 'w']
        self.patients['p10'].preferences = ['k11', 'k1', 'k4', 'k5', 'k6', 'k7','k2','w']
        self.patients['p11'].preferences = ['k3', 'k6', 'k5', 'k11']
        self.patients['p12'].preferences = ['k11', 'k3', 'k5', 'k9', 'k8', 'k10', 'k12']
        self._log("Paper example loaded successfully.")
        self.list_all_pairs()

# ########################################################################
# User Input and Main Menu
# ########################################################################
def _get_hla_from_user(person_type, locus):
    """Prompts the user to select two unique HLA types from a list for a given locus."""
    options = HLA_UNIVERSE[locus]
    print(f"\n--- Select 2 HLA-{locus} types for {person_type} ---")
    for i, hla in enumerate(options):
        print(f" {i+1:2}: {hla:3}", end='\n' if (i+1) % 8 == 0 else '\t')
    print()
    while True:
        input_str = input(f"Enter 2 numbers for HLA-{locus} (e.g., 1, 5): ")
        try:
            # 1. Attempt to convert input to a list of integers
            choices_str = [c.strip() for c in input_str.split(',')]
            # 2. Check for the correct number of inputs
            if len(choices_str) != 2:
                print("Error: Please enter exactly two numbers.")
                continue
            choices = [int(c) for c in choices_str]
            # 3. Check for duplicate inputs
            if len(set(choices)) != 2:
                print("Error: Please enter two unique (non-repeating) numbers.")
                continue
            # 4. Check if all numbers are within the valid range
            if not all(1 <= c <= len(options) for c in choices):
                print(f"Error: Numbers must be between 1 and {len(options)}.")
                continue
            # If all checks pass, return the selected HLA types
            return [options[c - 1] for c in choices]

        except ValueError:
            print("Error: Invalid input. Please use numbers separated by a comma.")

def _prompt_for_full_hla(person_type):
    return {locus: _get_hla_from_user(person_type, locus) for locus in ['A', 'B', 'DR']}

def _get_unacceptable_antigens_from_user():
    unacceptable = set()
    print("\n--- Enter Patient's Unacceptable Antigens (optional) ---")
    all_hlas = [hla for locus in sorted(HLA_UNIVERSE.keys()) for hla in HLA_UNIVERSE[locus]]
    for i, hla in enumerate(all_hlas): print(f" {i+1:2}: {hla:3}", end='\n' if (i+1)%8==0 else '\t')
    print()
    
    choice_str = input("Enter numbers separated by commas (or leave blank): ")
    if choice_str.strip():
        try:
            choices = {int(c.strip()) for c in choice_str.split(',')}
            for c in choices:
                if 1 <= c <= len(all_hlas):
                    unacceptable.add(all_hlas[c-1])
        except ValueError:
            print("Invalid input for unacceptable antigens, skipping.")
    return unacceptable

def main():
    ke = KidneyExchange()
    while True:
        print("\n===== Kidney Exchange TTCC Simulator Menu =====")
        print("1. Add Patient-Donor Pair")
        print("2. Generate All Preferences (incl. Crossmatch)")
        print("3. List All Pairs")
        print("4. Run TTCC Algorithm")
        print("5. Save Current State")
        print("6. Load State")
        print("7. Reset System")
        print("8. Load Paper's 12-Pair Example")
        print("0. Exit")
        choice = input("Enter your choice: ")
        
        if choice == '1':
            print("\n--- Adding New Patient-Donor Pair ---")
            try:
                patient_data = {
                    'blood_type': input("Patient blood type (O,A,B,AB): ").upper(),
                    'age': int(input("Patient age: ")),
                    'pra': int(input("Patient PRA (0-100): ")),
                    'hla': _prompt_for_full_hla("Patient"),
                    'unacceptable_antigens': _get_unacceptable_antigens_from_user()
                }
                donor_data = {
                    'blood_type': input("Donor blood type (O,A,B,AB): ").upper(),
                    'age': int(input("Donor age: ")),
                    'hla': _prompt_for_full_hla("Donor")
                }
                wants_w = input("Patient willing to consider waitlist? (y/n): ").lower() == 'y'
                ke.add_pair(patient_data, donor_data, wants_w)
            except ValueError:
                print("Invalid numerical input. Please try again.")

        elif choice == '2':
            ke.generate_all_preferences()
        
        elif choice == '3':
            ke.list_all_pairs()
            
        elif choice == '4':
            if not ke.patients: print("No pairs in system."); continue
            if not any(p.preferences for p in ke.patients.values()):
                print("Please generate preferences first (Option 2)."); continue

            try:
                max_c_str = input("Enter max cycle length (pairs) [Enter for no limit]: ")
                max_c = int(max_c_str) if max_c_str else 999
                max_ch_str = input("Enter max chain length (pairs) [Enter for no limit]: ")
                max_ch = int(max_ch_str) if max_ch_str else 999
            except ValueError:
                print("Invalid number. Using no limit."); max_c, max_ch = 999, 999

            print("\nSelect Chain Selection Rule:")
            print(" a: Min length, remove          | d: Priority, remove")
            print(" b: Max length, remove          | e: Priority, keep")
            print(" c: Max length, keep            | f: Hybrid O-type")
            print(" g: Best-Value (score-based)")
            rule = input("Enter rule (a-g): ").lower()
            if rule in ['a', 'b', 'c', 'd', 'e', 'f', 'g']:
                ke.run_ttcc(chain_rule=rule, max_cycle_len=max_c, max_chain_len=max_ch)
            else:
                print("Invalid rule.")

        elif choice == '5':
            ke.save_state(input("Enter filename to save (e.g., state.json): "))
        elif choice == '6':
            ke.load_state(input("Enter filename to load (e.g., state.json): "))
        elif choice == '7':
            ke.reset()
        elif choice == '8':
            ke._load_paper_example()
        elif choice == '0':
            print("Exiting simulator."); break
        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()