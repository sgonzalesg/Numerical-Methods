# Numerical Methods: Scientific Computing & Skill Development

Welcome to my Numerical Methods repository! This space is dedicated to deepening my understanding, implementation, and optimization of essential numerical algorithms used to solve complex mathematical and scientific problems.

## Project Purpose

The primary objective of this repository is to **strengthen my scientific computing and software engineering skills**. By translating theoretical mathematical frameworks into functional, modular code, I aim to master algorithmic efficiency, error analysis, and the convergence behavior of numerical techniques.

---

## Repository Structure

The project is organized by mathematical domains. Each module contains the core algorithmic logic along with test scripts and demonstration notebooks:

*   **`equations/` (Roots of Non-Linear Equations):** Bracketing methods (Bisection) and open methods (Newton-Raphson, Secant) to solve $f(x) = 0$.
*   **`linear_algebra/` (Systems of Equations & Matrices):** Direct solvers (LU Decomposition, Gauss-Jordan) and iterative methods (Jacobi, Gauss-Seidel) for $Ax = b$.
*   **`integration/` (Numerical Integration):** Newton-Cotes formulas, including the Trapezoidal Rule and Simpson's rules (1/3 and 3/8).
*   **`ode/` (Ordinary Differential Equations):** Initial value problem solvers using single-step methods like Euler and 4th-Order Runge-Kutta (RK4).

---

## Tech Stack & Libraries

To maintain industrial and academic standards in scientific computing, this project utilizes:

*   **Python 3.x:** The core programming language chosen for its readability and robust scientific ecosystem.
*   **NumPy:** For high-performance vector and matrix operations, ensuring memory-efficient linear algebra.
*   **Matplotlib:** To visualize method convergence, error decay, and geometric approximations.
*   **Jupyter Notebooks:** Used as an interactive sandbox environment for experimentation and detailed data visualization.

---

## Getting Started

### Prerequisites

Ensure you have Python 3.8+ installed on your system along with a virtual environment manager (`venv` or `conda`).

### Installation & Setup

1. Clone the repository to your local machine:
```bash
   git clone [https://github.com/YOUR_USERNAME/numerical-methods.git](https://github.com/YOUR_USERNAME/numerical-methods.git)
   cd numerical-methods
```

2. Create and activate a virtual environment:
```bash
   python -m venv venv
   # On Windows:
   .\venv\Scripts\activate
   # On Linux/macOS:
   source venv/bin/activate
```

3. Install the required dependencies:
```bash
   pip install -r requirements.txt
```
