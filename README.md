<p align="center">
	<img src="./asset/logo.png" alt="FMAD logo" width="600" height="200">
</p>

<h1 align="center">FMAD</h1>

<p align="center">
	<strong>Calculate derivatives using Forward Mode Automatic Differentiation!</strong>
</p>

> [!NOTE]
> The `normal` method for calculating the derivative of an equation only uses the power rule currently, and not a rule-based methodology that would include the sum, product, quotient, and chain rules.

## 🚀 Overview

Welcome to **FMAD**! This program allows you to calculate the derivatives of a function using Forward Mode Automatic Differentiation. The program is designed to be user-friendly and easy to use. The program is written in C++.

## 🎨 Features

- **Multiple Methods:** The program uses the `normal` method, limit definition, and forward mode automatic differentiation to calculate the derivative of a function.
- **Key Values:** The program outputs key values (important points on the graph that would show up in a data table) for the given equation and the derivative equation.
- **Comparing Methods:** The program displays from fastest to slowest, the times it takes for each method to compute the derivative of the function, with color interpolation from green to red showing how fast the method took in comparison to the average time of all of the methods.

## 🛠️ Installation

To get started with the program, follow the steps below:

1. **Clone the Repository**
```sh
git clone https://github.com/321BadgerCode/fmad.git
cd ./fmad/
```

2. **Compile the Program**
```sh
g++ ./main.cpp -o fmad
```

## 📈 Usage

To use the program, there is only **one** step!

1. **Run the program**
```sh
./fmad <equation>
```

<details>

<summary>Examples</summary>

```sh
./fmad "x^2"
```

```sh
./fmad "x^3 + 2x^2 + 3x + 4"
```

</details>

## 📝 Documentation

```sh
man ./fmad.1
```

## 📜 License

[LICENSE](./LICENSE)