# Exam Preparation Guide - Constraint Programming

Use this checklist and topic summary to prepare for the oral exam.

## Core Concepts to Master
### 1. Modeling in MiniZinc
- [ ] **Variable Domains:** Defining appropriate domains (int, bool, sets).
- [ ] **Global Constraints:** Understanding when to use `alldifferent`, `cumulative`, `circuit`, `table`.
- [ ] **Reification:** Using boolean variables to represent the satisfaction of constraints (crucial for complex logic).
- [ ] **Redundant Constraints:** Adding constraints that don't change the solution set but help propagation.

### 2. Consistency & Propagation
- [ ] **Arc Consistency (AC):** AC1, AC3. How they work to prune domains.
- [ ] **Bounds Consistency:** Used for numerical constraints.
- [ ] **Propagator Filtering:** The mechanism behind `alldifferent` (Regin's theorem, bipartite graphs).

### 3. Search Strategies
- [ ] **Search Tree:** Understanding nodes, branches, and depth.
- [ ] **Variable Selection:** `first_fail`, `most_constrained`, etc.
- [ ] **Domain Selection:** `indomain_min`, `indomain_random`, etc.
- [ ] **Optimization:** Branch and Bound for minimization/maximization.

### 4. Symmetry Breaking
- [ ] **Identifying Symmetries:** Rotation, reflection, variable swapping.
- [ ] **Symmetry Breaking Constraints:** Lexicographic order constraints to remove redundant branches.

## Exam To-Do List
- [ ] **Finalize Project:** Ensure your report is complete and you can explain every line of your code.
- [ ] **Theory Review:** Be ready to explain AC3 or Regin's theorem on a whiteboard/digital screen.
- [ ] **Problem Portfolio:** Review the models for N-Queens, Sudoku, Timetabling, and TSP.
- [ ] **(Optional) Paper Review:** If you chose to review a paper (e.g., Nurse Scheduling), prepare a 5-10 minute summary.

## Key Theoretical Theorems
- **Berge's Theorem:** Relating matchings to alternating cycles in bipartite graphs.
- **Regin's Theorem:** For filtering the `alldifferent` constraint.
- **NP-Completeness:** Why CSP is generally NP-complete and how we handle it in practice.

---
*Good luck with your exam!*
