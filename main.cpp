#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <cctype>
#include <ctime>
#include <cmath>
#include <functional>

// Structure for polynomial terms
struct Term {
	double coefficient;
	int exponent;
};

// Dual number struct
struct Dual {
	double real;
	double dual;

	Dual(double r, double d) : real(r), dual(d) {}

	// Overload addition
	Dual operator+(const Dual &other) const {
		return Dual(real + other.real, dual + other.dual);
	}

	// Overload multiplication
	Dual operator*(const Dual &other) const {
		return Dual(real * other.real, real * other.dual + dual * other.real);
	}

	// Overload exponentiation for polynomials (x^n)
	Dual pow(int exponent) const {
		if (exponent == 0) return Dual(1, 0);
		if (exponent == 1) return *this;
		Dual result = *this;
		for (int i = 1; i < exponent; i++) {
			result = result * (*this);
		}
		return result;
	}
};

// Structure for rgb
struct RGB {
	int r;
	int g;
	int b;

	RGB(int red, int green, int blue) : r(red), g(green), b(blue) {}
	RGB(std::vector<int> colors) {
		r = colors[0];
		g = colors[1];
		b = colors[2];
	}

	RGB operator +(const RGB &other) const {
		return RGB(r + other.r, g + other.g, b + other.b);
	}
	RGB operator -(const RGB &other) const {
		return RGB(r - other.r, g - other.g, b - other.b);
	}
	RGB operator *(float scalar) const {
		return RGB(r * scalar, g * scalar, b * scalar);
	}

	RGB interpolate(RGB end, float percent) {
		return (*this + (end - *this) * percent);
	}
};

// Replace all function
void replaceAll(std::string &str, const std::string &from, const std::string &to) {
	size_t start_pos = 0;
	while ((start_pos = str.find(from, start_pos)) != std::string::npos) {
		str.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
}

// Function to parse a polynomial expression
// TODO: Add support for more complex expressions (e.g. sin, cos, etc.) and parentheses by returning a tree structure from tokenized input
std::vector<Term> parsePolynomial(const std::string &expr) {
	std::string preprocessed_expr = expr;
	replaceAll(preprocessed_expr, " ", "");
	std::vector<Term> terms;
	std::istringstream iss(preprocessed_expr);
	std::string token;
	size_t i = 0;

	while (i < preprocessed_expr.length()) {
		double coefficient = 1.0;
		int exponent = 0;
		bool negative = false;

		// Handle sign
		if (preprocessed_expr[i] == '-') {
			negative = true;
			i++;
		} else if (preprocessed_expr[i] == '+') {
			i++;
		}

		// Read coefficient
		size_t start = i;
		while (i < preprocessed_expr.length() && (isdigit(preprocessed_expr[i]) || preprocessed_expr[i] == '.')) i++;
		if (start != i) {  // Number found
			coefficient = std::stod(preprocessed_expr.substr(start, i - start));
		}

		if (i < preprocessed_expr.length() && preprocessed_expr[i] == 'x') { // Found 'x'
			i++; // Move past 'x'
			exponent = 1; // Default exponent is 1 if no caret found

			if (i < preprocessed_expr.length() && preprocessed_expr[i] == '^') { // Found '^'
				i++; // Move past '^'
				start = i;
				while (i < preprocessed_expr.length() && isdigit(preprocessed_expr[i])) i++;
				exponent = std::stoi(preprocessed_expr.substr(start, i - start));
			}
		}

		if (negative) coefficient = -coefficient;
		terms.push_back({coefficient, exponent});
	}

	return terms;
}

// Compute derivative using power rule
std::vector<Term> differentiate(const std::vector<Term> &terms) {
	// TODO: Implement rule based methodology, including:
		// Sum rule: d/dx (f(x) + g(x)) = f'(x) + g'(x)
		// Product rule: d/dx (f(x) * g(x)) = f'(x) * g(x) + f(x) * g'(x)
		// Quotient rule: d/dx (f(x) / g(x)) = (f'(x) * g(x) - f(x) * g'(x)) / g(x)^2
		// Chain rule: d/dx (f(g(x))) = f'(g(x)) * g'(x)
	std::vector<Term> derivativeTerms;

	for (const auto &term : terms) {
		if (term.exponent == 0) continue; // Skip constants (derivative is 0)

		double coeff = term.coefficient * term.exponent;
		int exponent = term.exponent - 1;

		derivativeTerms.push_back({coeff, exponent});
	}

	return derivativeTerms;
}

// Compute derivative using dual numbers
std::vector<Term> differentiateWithDuals(const std::vector<Term> &terms) {
	std::vector<Term> derivativeTerms;
	Dual x(1.0, 1.0); // Set x = (1 + ε) for AD

	for (const auto &term : terms) {
		if (term.exponent == 0) continue; // Skip constants (derivative is 0)

		Dual coeff(term.coefficient, 0);
		Dual termResult = coeff * x.pow(term.exponent);

		derivativeTerms.push_back({termResult.dual, term.exponent - 1});
	}

	return derivativeTerms;
}

// Compute derivative using limits
std::vector<Term> differentiateWithLimits(const std::vector<Term> &terms, double h = 1e-10) {
	std::vector<Term> derivativeTerms;

	for (const auto &term : terms) {
		if (term.exponent == 0) continue; // Skip constants (derivative is 0)

		double coeff = term.coefficient;
		int exponent = term.exponent;

		double derivative = (coeff * std::pow(1 + h, exponent) - coeff) / h;
		derivativeTerms.push_back({derivative, exponent - 1});
	}

	return derivativeTerms;
}

// Convert terms back to string representation
std::string termsToString(const std::vector<Term> &terms) {
	if (terms.empty()) return "0"; // Edge case: constant function

	std::ostringstream oss;
	bool first = true;

	for (const auto &term : terms) {
		if (!first) {
			if (term.coefficient > 0) oss << " + ";
			else oss << " - ";
		} else if (term.coefficient < 0) {
			oss << "-";
		}

		double absCoeff = std::abs(term.coefficient);
		if (!(absCoeff == 1.0 && term.exponent > 0)) oss << absCoeff; // Omit '1' for terms like 'x^n'

		if (term.exponent > 0) oss << "x";
		if (term.exponent > 1) oss << "^" << term.exponent;

		first = false;
	}
	return oss.str();
}

// Function to take in string function and return the c++ function for it
std::function<double(double)> getFunction(const std::string &expr) {
	std::vector<Term> terms = parsePolynomial(expr);

	return [terms](double x) {
		double result = 0;
		for (const auto &term : terms) {
			result += term.coefficient * std::pow(x, term.exponent);
		}
		return result;
	};
}

// Function to get key values of function (points that go in a table to best express the equation) by using critical points, inflection points, zeroes, and endpoints
std::vector<std::pair<double, double>> getKeyValues(const std::vector<Term> &terms) {
	std::vector<std::pair<double, double>> keyValues;
	std::function<double(double)> f = getFunction(termsToString(terms));
	std::function<double(double)> fPrime = getFunction(termsToString(differentiate(terms)));

	// Critical points
	std::vector<double> criticalPoints;
	for (int i = 0; i < terms.size(); i++) {
		if (terms[i].exponent == 1) {
			criticalPoints.push_back(-terms[i].coefficient / terms[i + 1].coefficient);
		}
	}
	for (double x : criticalPoints) {
		keyValues.push_back({x, f(x)});
	}

	// Inflection points
	std::vector<double> inflectionPoints;
	std::vector<Term> secondDerivative = differentiate(differentiate(terms));
	std::function<double(double)> fDoublePrime = getFunction(termsToString(secondDerivative));
	for (int i = 0; i < terms.size(); i++) {
		if (terms[i].exponent == 2) {
			inflectionPoints.push_back(-terms[i].coefficient / terms[i + 1].coefficient);
		}
	}
	for (double x : inflectionPoints) {
		keyValues.push_back({x, f(x)});
	}

	// Zeroes
	std::vector<double> zeroes;
	for (int i = 0; i < terms.size(); i++) {
		if (terms[i].coefficient == 0) {
			zeroes.push_back(0);
		}
	}
	for (double x : zeroes) {
		keyValues.push_back({x, f(x)});
	}

	// Endpoints
	std::vector<double> endpoints = {0, 1, -1};
	for (double x : endpoints) {
		keyValues.push_back({x, f(x)});
	}

	// Delete duplicates
	for (int i = 0; i < keyValues.size(); i++) {
		for (int j = i + 1; j < keyValues.size(); j++) {
			if (keyValues[i].first == keyValues[j].first) {
				keyValues.erase(keyValues.begin() + j);
			}
		}
	}

	return keyValues;
}

// Function to sort numbers from least to greatest
std::vector<double> sortNumbers(std::vector<double> numbers) {
	for (int i = 0; i < numbers.size(); i++) {
		for (int j = i + 1; j < numbers.size(); j++) {
			if (numbers[i] > numbers[j]) {
				double temp = numbers[i];
				numbers[i] = numbers[j];
				numbers[j] = temp;
			}
		}
	}
	return numbers;
}

// Function to return average of numbers
double average(std::vector<double> numbers) {
	double sum = 0;
	for (int i = 0; i < numbers.size() - 1; i++) {
		sum += numbers[i];
	}
	return sum / numbers.size();
}

// Function to print a given string with unicode color of rgb
std::string getColorTxt(std::string text, RGB rgb) {
	std::string color = "\033[38;2;" + std::to_string(rgb.r) + ";" + std::to_string(rgb.g) + ";" + std::to_string(rgb.b) + "m";
	return color + text + "\033[0m";
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		std::cerr << "Usage: " << argv[0] << " \"expression\"" << std::endl;
		return 1;
	}

	std::string equation = argv[1];

	std::vector<Term> terms = parsePolynomial(equation);

	clock_t start = clock();
	std::vector<Term> derivative = differentiate(terms);
	double end = (clock() - start) / (double)CLOCKS_PER_SEC;

	start = clock();
	std::vector<Term> derivativeWithDuals = differentiateWithDuals(terms);
	double end2 = (clock() - start) / (double)CLOCKS_PER_SEC;

	start = clock();
	std::vector<Term> derivativeWithLimits = differentiateWithLimits(terms);
	double end3 = (clock() - start) / (double)CLOCKS_PER_SEC;

	std::cout << "f(x) = " << equation << std::endl;
	std::vector<std::pair<double, double>> keyValues = getKeyValues(terms);
	std::cout << "Key values:" << std::endl;
	for (const auto &pair : keyValues) {
		std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
	}
	std::cout << std::endl;
	std::cout << "f'(x) = " << termsToString(derivativeWithDuals) << std::endl;
	keyValues = getKeyValues(derivativeWithDuals);
	std::cout << "Key values:" << std::endl;
	for (const auto &pair : keyValues) {
		std::cout << "(" << pair.first << ", " << pair.second << ")" << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Time taken for each method:" << std::endl;
	std::vector<double> times = {end, end2, end3};
	std::vector<double> sortedTimes = sortNumbers(times);
	double avg = average(sortedTimes);
	std::vector<bool> used = {false, false, false};
	for (int i = 0; i < sortedTimes.size(); i++) {
		double percentDiff = 100 - ((avg / sortedTimes[i]) * 100);
		RGB color = RGB(0, 255, 0).interpolate(RGB(255, 0, 0), percentDiff / 100);
		if (sortedTimes[i] == end && !used[0]) {
			std::cout << "Normal derivative: " << getColorTxt(std::to_string(end) + "s", color) << std::endl;
			used[0] = true;
		} else if (sortedTimes[i] == end2 && !used[1]) {
			std::cout << "Derivative with duals: " << getColorTxt(std::to_string(end2) + "s", color) << std::endl;
			used[1] = true;
		} else if (sortedTimes[i] == end3 && !used[2]) {
			std::cout << "Derivative with limits: " << getColorTxt(std::to_string(end3) + "s", color) << std::endl;
			used[2] = true;
		}
	}

	return 0;
}