#pragma once
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "constants.h"

// Class that performs linear algebraic opertations on 3-dimensional vectors
class Vector3 {
private:
	double m_x {};
	double m_y {};
	double m_z {};
	
public:
	constexpr Vector3(double x, double y, double z) // defines new variable type Vector3, which takes in arguments (x, y, z)
		: m_x{ x }, m_y{ y }, m_z{ z }{}
	constexpr Vector3() : m_x{ 0.0 }, m_y{ 0.0 }, m_z{ 0.0 } {} // sets the default to (0,0,0)

	constexpr Vector3 operator*(double scalar) const
	{
		return Vector3(m_x * scalar, m_y * scalar, m_z * scalar); // multiplies a scalar into a Vector3
	}

	constexpr Vector3 operator+(const Vector3& u) const // adds 2 Vector3's
	{
		return Vector3(m_x + u.m_x, m_y + u.m_y, m_z + u.m_z);
	}

	constexpr Vector3 operator-(const Vector3& u) const // subtracts 2 Vector3's
	{
		return Vector3(m_x - u.m_x, m_y - u.m_y, m_z - u.m_z);
	}

	constexpr Vector3 mult(const Vector3& u) const // multiplies 2 Vector3's
	{
		return Vector3(m_x * u.m_x, m_y * u.m_y, m_z * u.m_z);
	}

	double magnitude() const // finds the magnitude of a Vector3
	{
		return {std::sqrt( (m_x * m_x) + (m_y * m_y) + (m_z * m_z) )}; // returns a double
	}

	Vector3 normalize() const // normalize inputted Vector3
	{
		double mag { magnitude() }; // determine magnitude of input
		if (mag == 0.0)
			return { 0.0, 0.0, 0.0 }; // in case of zero magnitude, to avoid undefined variables
		return { m_x / mag, m_y / mag, m_z / mag };
	}

	constexpr double dotP(const Vector3& u) const // dot-product operation using an additional Vector3 type
	{
		return (m_x * u.m_x) + (m_y * u.m_y) + (m_z * u.m_z);
	}

	constexpr Vector3 crossP(const Vector3& u) const // cross-product operation using an additional Vector3 type
	{
		double cx = (m_y * u.m_z) - (m_z * u.m_y);
		double cy = (m_z * u.m_x) - (m_x * u.m_z);
		double cz = (m_x * u.m_y) - (m_y * u.m_x);

		return Vector3(cx, cy, cz);
	}

	// extracts individual components of Vector3 
	double getX() const { return m_x; }
	double getY() const { return m_y; }
	double getZ() const { return m_z; }

	void print() const
	{
		std::cout << "(" << m_x << ", " << m_y << ", " << m_z << ")\n";
	}

	friend std::ostream& operator<<(std::ostream& out, const Vector3& v);

};

inline std::ostream& operator<<(std::ostream& out, const Vector3& v)
{
	out << "(" << v.m_x << ", " << v.m_y << ", " << v.m_z << ")";
	return out;
}

inline constexpr Vector3 operator*(double scalar, const Vector3& v)
{
	return v * scalar; // flips order of scalar multiplication
}

// Class that performs linear algebraic operations on 3x3 matrices, which are used for rotation transformations
class Matrix3x3 { 
private:
	// Each row of the matrix is represented as a Vector3
	Vector3 m_row0{};
	Vector3 m_row1{};
	Vector3 m_row2{};
public:
	// Constructor that takes in three Vector3s as rows of the matrix
	constexpr Matrix3x3(const Vector3& row0, const Vector3& row1, const Vector3& row2) 
		: m_row0{ row0 }, m_row1{ row1 }, m_row2{ row2 } {}
	constexpr Matrix3x3() : m_row0{ 1.0, 0.0, 0.0 }, m_row1{ 0.0, 1.0, 0.0 }, m_row2{ 0.0, 0.0, 1.0 } {} // default to identity matrix

	// Finds the transpose of the matrix
	constexpr Matrix3x3 trans() const
	{
		Vector3 newRow0{ m_row0.getX(), m_row1.getX(), m_row2.getX() };
		Vector3 newRow1{ m_row0.getY(), m_row1.getY(), m_row2.getY() };
		Vector3 newRow2{ m_row0.getZ(), m_row1.getZ(), m_row2.getZ() };
		return Matrix3x3(newRow0, newRow1, newRow2);
	}

	// Multiplies the matrix by a scalar, resulting in a new Matrix3x3
	constexpr Matrix3x3 operator*(double scalar) const
	{
		return Matrix3x3(m_row0 * scalar, m_row1 * scalar, m_row2 * scalar);
	}

	constexpr Matrix3x3 operator+(const Matrix3x3& u) const // adds two matrices together
	{
		return Matrix3x3(m_row0 + u.m_row0, m_row1 + u.m_row1, m_row2 + u.m_row2);
	}

	constexpr Matrix3x3 operator-(const Matrix3x3& u) const // subtracts two matrices together
	{
		return Matrix3x3(m_row0 - u.m_row0, m_row1 - u.m_row1, m_row2 - u.m_row2);
	}

	// Multiplies the matrix by a Vector3, resulting in a new Vector3
	constexpr Vector3 operator*(const Vector3& v) const
	{
		double x = m_row0.dotP(v);
		double y = m_row1.dotP(v);
		double z = m_row2.dotP(v);
		return Vector3(x, y, z);
	}

	// Multiplies the matrix by another Matrix3x3, resulting in a new Matrix3x3
	constexpr Matrix3x3 operator*(const Matrix3x3& u) const
	{
		Matrix3x3 uT{ u.trans() }; // transpose the second matrix to facilitate multiplication
		return Matrix3x3(
			(*this) * uT.m_row0, // multiply this matrix by each row of the other transposed matrix
			(*this) * uT.m_row1,
			(*this) * uT.m_row2
		).trans(); // transpose back to get correct rows
	}

	constexpr double det() const
	{
		double a = m_row0.getX(), b = m_row0.getY(), c = m_row0.getZ();
		double d = m_row1.getX(), e = m_row1.getY(), f = m_row1.getZ();
		double g = m_row2.getX(), h = m_row2.getY(), i = m_row2.getZ();
		double determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
		return determinant;
	}

	constexpr Matrix3x3 inv() const
	{
		double determinant = det();
		if (determinant == 0.0)
			throw std::runtime_error("Matrix is singular and cannot be inverted.");
		double a = m_row0.getX(), b = m_row0.getY(), c = m_row0.getZ();
		double d = m_row1.getX(), e = m_row1.getY(), f = m_row1.getZ();
		double g = m_row2.getX(), h = m_row2.getY(), i = m_row2.getZ();
		Matrix3x3 adjugate(
			Vector3(e * i - f * h, c * h - b * i, b * f - c * e),
			Vector3(f * g - d * i, a * i - c * g, c * d - a * f),
			Vector3(d * h - e * g, b * g - a * h, a * e - b * d)
		);
		return adjugate * (1.0 / determinant);
	}
};


inline constexpr Matrix3x3 operator*(double scalar, const Matrix3x3& v)
{
	return v * scalar; // flips order of scalar multiplication
}

// Class that performs linear algebraic operations on quaternions, which are used for rotation transformations in 3D space. 
class Quaternion
{
private: // Quaternions are represented as (w, x, y, z), where w is the scalar part and (x, y, z) is the vector part.
	double m_w{};
	double m_x{};
	double m_y{};
	double m_z{};
public:
	constexpr Quaternion(double w, double x, double y, double z)
		: m_w{ w }, m_x{ x }, m_y{ y }, m_z{ z } {}
	constexpr Quaternion() : m_w{ 1.0 }, m_x{ 0.0 }, m_y{ 0.0 }, m_z{ 0.0 } {} // default to identity quaternion

	double getW() const { return m_w; }
	double getX() const { return m_x; }
	double getY() const { return m_y; }
	double getZ() const { return m_z; }

	constexpr Quaternion operator*(double scalar) const
	{
		return Quaternion(m_w * scalar, m_x * scalar, m_y * scalar, m_z * scalar); // multiplies a scalar into a quaternion
	}

	constexpr Quaternion operator+(const Quaternion& q) const // adds two quaternions together
	{
		return Quaternion(m_w + q.m_w, m_x + q.m_x, m_y + q.m_y, m_z + q.m_z);
	}

	constexpr Quaternion operator*(const Quaternion& q) const // multiplies two quaternions together, which is used to combine rotations
	{
		double w = m_w * q.m_w - m_x * q.m_x - m_y * q.m_y - m_z * q.m_z;
		double x = m_w * q.m_x + m_x * q.m_w + m_y * q.m_z - m_z * q.m_y;
		double y = m_w * q.m_y - m_x * q.m_z + m_y * q.m_w + m_z * q.m_x;
		double z = m_w * q.m_z + m_x * q.m_y - m_y * q.m_x + m_z * q.m_w;
		return Quaternion(w, x, y, z);
	}

	constexpr Quaternion conj() const // returns the conjugate of the quaternion, which is used in rotation operations
	{
		return Quaternion(m_w, -m_x, -m_y, -m_z);
	}

	double magnitude() const // finds the magnitude of the quaternion, which is used for normalization
	{
		return std::sqrt(m_w * m_w + m_x * m_x + m_y * m_y + m_z * m_z);
	}

	Quaternion normalize() const // normalizes the quaternion to ensure it represents a valid rotation (unit quaternion)
	{
		double mag = magnitude();
		if (mag == 0.0)
			return Quaternion(1.0, 0.0, 0.0, 0.0); // in case of zero magnitude, to avoid undefined variables
		return Quaternion(m_w / mag, m_x / mag, m_y / mag, m_z / mag);
	}

	constexpr Matrix3x3 toMatrix() const // converts the quaternion to a corresponding rotation matrix, which can be used to rotate vectors in 3D space
	{
		double w = m_w, x = m_x, y = m_y, z = m_z;
		return Matrix3x3(
			Vector3(1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)),
			Vector3(2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)),
			Vector3(2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y))
		);
	}

	// Creates a quaternion from an axis of rotation and an angle in degrees, which can be used to represent a rotation in 3D space
	constexpr static Quaternion fromAxisAngle(const Vector3& axis, double angle)
	{
		Vector3 normAxis{ axis.normalize() };
		double rad = angle * Constants::pi / 180.0; // convert to radians
		double halfAngle = rad / 2.0;
		double s = std::sin(halfAngle);
		return Quaternion(std::cos(halfAngle), normAxis.getX() * s, normAxis.getY() * s, normAxis.getZ() * s).normalize();
	}

	// Rotates a Vector3 using the quaternion, which is done by treating the vector as a quaternion with a scalar part of 0 and applying the rotation formula: q * v * q^-1
	Vector3 rotate(const Vector3& v) const
	{
		Quaternion p(0, v.getX(), v.getY(), v.getZ());
		Quaternion result = (*this) * p * this->conj();
		return Vector3(result.m_x, result.m_y, result.m_z);
	}
};
inline constexpr Quaternion operator*(double scalar, const Quaternion& v)
{
	return v * scalar; // flips order of scalar multiplication
}


// Class that performs linear algebraic operations on 6-dimensional vectors, which are used for state vectors in the projectile motion simulation (position and velocity in 3D space)
class Vector6 {
private:
	double m_x1{};
	double m_x2{};
	double m_x3{};
	double m_x4{};
	double m_x5{};
	double m_x6{};
public:
	constexpr Vector6(double x1, double x2, double x3, double x4, double x5, double x6)
		: m_x1{ x1 }, m_x2{ x2 }, m_x3{ x3 }, m_x4{ x4 }, m_x5{ x5 }, m_x6{ x6 } {
	}
	constexpr Vector6() : m_x1{ 0.0 }, m_x2{ 0.0 }, m_x3{ 0.0 }, m_x4{ 0.0 }, m_x5{ 0.0 }, m_x6{ 0.0 } {} // default to zero vector

	constexpr Vector6 operator*(double scalar) const
	{
		return Vector6(m_x1 * scalar, m_x2 * scalar, m_x3 * scalar, m_x4 * scalar, m_x5 * scalar, m_x6 * scalar); // multiplies a scalar into a Vector6
	}

	constexpr Vector6 operator+(const Vector6& u) const // adds 2 Vector6's
	{
		return Vector6(m_x1 + u.m_x1, m_x2 + u.m_x2, m_x3 + u.m_x3, m_x4 + u.m_x4, m_x5 + u.m_x5, m_x6 + u.m_x6);
	}

	constexpr Vector6 operator-(const Vector6& u) const // subtracts 2 Vector6's
	{
		return Vector6(m_x1 - u.m_x1, m_x2 - u.m_x2, m_x3 - u.m_x3, m_x4 - u.m_x4, m_x5 - u.m_x5, m_x6 - u.m_x6);
	}

	constexpr double dotP(const Vector6& u) const // dot-product operation using an additional Vector6 type
	{
		return (m_x1 * u.m_x1) + (m_x2 * u.m_x2) + (m_x3 * u.m_x3) + (m_x4 * u.m_x4) + (m_x5 * u.m_x5) + (m_x6 * u.m_x6);
	}

	double magnitude() const // finds the magnitude of a Vector6
	{
		return { std::sqrt((m_x1 * m_x1) + (m_x2 * m_x2) + (m_x3 * m_x3) + (m_x4 * m_x4) + (m_x5 * m_x5) + (m_x6 * m_x6)) }; // returns a double
	}

	double getX1() const { return m_x1; }
	double getX2() const { return m_x2; }
	double getX3() const { return m_x3; }
	double getX4() const { return m_x4; }
	double getX5() const { return m_x5; }
	double getX6() const { return m_x6; }

	double get(int i) const // returns component i
	{
		switch (i) {
		case 0: return m_x1;
		case 1: return m_x2;
		case 2: return m_x3;
		case 3: return m_x4;
		case 4: return m_x5;
		case 5: return m_x6;
		default: throw std::out_of_range("Index out of range for Vector6");
		}
	}

	Vector6 normalize() const // normalize inputted Vector6
	{
		double mag{ magnitude() }; // determine magnitude of input
		if (mag == 0.0)
			return { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }; // in case of zero magnitude, to avoid undefined variables
		return { m_x1 / mag, m_x2 / mag, m_x3 / mag, m_x4 / mag, m_x5 / mag, m_x6 / mag };
	}

	friend std::ostream& operator<<(std::ostream& out, const Vector6& v); 
};

inline std::ostream& operator<<(std::ostream& out, const Vector6& v)
{
	out << "(" << v.m_x1 << ", " << v.m_x2 << ", " << v.m_x3 << ", " << v.m_x4 << ", " << v.m_x5 << ", " << v.m_x6 << ")";
	return out;
}

// Class that performs linear algebraic operations on 6x6 matrices, which are used for state transition matrices in the projectile motion simulation
class Matrix6x6 {
private:
	Vector6 m_row0{};
	Vector6 m_row1{};
	Vector6 m_row2{};
	Vector6 m_row3{};
	Vector6 m_row4{};
	Vector6 m_row5{};
public:
	constexpr Matrix6x6(const Vector6& row0, const Vector6& row1, const Vector6& row2, const Vector6& row3, const Vector6& row4, const Vector6& row5)
		: m_row0{ row0 }, m_row1{ row1 }, m_row2{ row2 }, m_row3{ row3 }, m_row4{ row4 }, m_row5{ row5 } {
	}
	constexpr Matrix6x6() : m_row0{ 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 }, m_row1{ 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 }, m_row2{ 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 }, m_row3{ 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 }, m_row4{ 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 }, m_row5{ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0 } {} // default to identity matrix

	constexpr Matrix6x6 operator*(double scalar) const
	{
		return Matrix6x6(m_row0 * scalar, m_row1 * scalar, m_row2 * scalar, m_row3 * scalar, m_row4 * scalar, m_row5 * scalar);
	}

	constexpr Matrix6x6 operator+(const Matrix6x6& u) const // adds two matrices together
	{
		return Matrix6x6(m_row0 + u.m_row0, m_row1 + u.m_row1, m_row2 + u.m_row2, m_row3 + u.m_row3, m_row4 + u.m_row4, m_row5 + u.m_row5);
	}

	constexpr Matrix6x6 operator-(const Matrix6x6& u) const // subtracts two matrices together
	{
		return Matrix6x6(m_row0 - u.m_row0, m_row1 - u.m_row1, m_row2 - u.m_row2, m_row3 - u.m_row3, m_row4 - u.m_row4, m_row5 - u.m_row5);
	}

	constexpr Vector6 operator*(const Vector6& v) const
	{
		double x1 = m_row0.dotP(v);
		double x2 = m_row1.dotP(v);
		double x3 = m_row2.dotP(v);
		double x4 = m_row3.dotP(v);
		double x5 = m_row4.dotP(v);
		double x6 = m_row5.dotP(v);
		return Vector6(x1, x2, x3, x4, x5, x6);
	}

	constexpr Matrix6x6 transpose() const
	{
		Vector6 newRow0{ m_row0.getX1(), m_row1.getX1(), m_row2.getX1(), m_row3.getX1(), m_row4.getX1(), m_row5.getX1() };
		Vector6 newRow1{ m_row0.getX2(), m_row1.getX2(), m_row2.getX2(), m_row3.getX2(), m_row4.getX2(), m_row5.getX2() };
		Vector6 newRow2{ m_row0.getX3(), m_row1.getX3(), m_row2.getX3(), m_row3.getX3(), m_row4.getX3(), m_row5.getX3() };
		Vector6 newRow3{ m_row0.getX4(), m_row1.getX4(), m_row2.getX4(), m_row3.getX4(), m_row4.getX4(), m_row5.getX4() };
		Vector6 newRow4{ m_row0.getX5(), m_row1.getX5(), m_row2.getX5(), m_row3.getX5(), m_row4.getX5(), m_row5.getX5() };
		Vector6 newRow5{ m_row0.getX6(), m_row1.getX6(), m_row2.getX6(), m_row3.getX6(), m_row4.getX6(), m_row5.getX6() };
		return Matrix6x6(newRow0, newRow1, newRow2, newRow3, newRow4, newRow5);
	}

	constexpr Matrix6x6 operator*(const Matrix6x6& u) const
	{
		Matrix6x6 uT{ u.transpose() }; // transpose the second matrix to facilitate multiplication
		return Matrix6x6(
			(*this) * uT.m_row0, // multiply this matrix by each row of the other transposed matrix
			(*this) * uT.m_row1,
			(*this) * uT.m_row2,
			(*this) * uT.m_row3,
			(*this) * uT.m_row4,
			(*this) * uT.m_row5
		).transpose(); // transpose back to get correct rows
	}
	
	Vector6 getRow(int i) const // returns row i as a Vector6
	{
		switch (i) {
		case 0: return m_row0;
		case 1: return m_row1;
		case 2: return m_row2;
		case 3: return m_row3;
		case 4: return m_row4;
		case 5: return m_row5;
		default: throw std::out_of_range("Index out of range for Matrix6x6");
		}
	}

	Matrix6x6 inv() const
	{
		// Done using Gaussian Elimination
		double aug[6][12]{}; // augmented matrix for Gaussian elimination

		// fill left 6 columns from this matrix
		// fill right 6 columns as identity
		for (int i{ 0 }; i < 6; ++i) {
			Vector6 row{ getRow(i) };
			for (int j{ 0 }; j < 6; ++j)
				aug[i][j] = row.get(j);      // left side -- this matrix
			aug[i][i + 6] = 1.0;             // right side -- identity
		}

		for (int k{ 0 }; k < 6; ++k)
		{
			// Find pivot for column k
			int maxRow = k;												// start with current row as pivot
			for (int i{ k + 1 }; i < 6; ++i) {							// check rows below for larger pivot
				if (std::abs(aug[i][k]) > std::abs(aug[maxRow][k]))		// if this row has a larger absolute value in column k, it becomes the new pivot
					maxRow = i;
			}

			if (std::abs(aug[maxRow][k]) < 1e-12)						// if the pivot is effectively zero, the matrix is singular and cannot be inverted
				throw std::runtime_error("Matrix6x6 is singular and cannot be inverted.");

			// Swap maximum row with current row (pivoting)
			for (int j{ 0 }; j < 12; ++j) {
				std::swap(aug[k][j], aug[maxRow][j]);					// swap the entire row in the augmented matrix to move the pivot into position
			}
			// Make all rows below this one 0 in current column
			for (int i{ k + 1 }; i < 6; ++i) {
				double factor = aug[i][k] / aug[k][k];					// calculate the factor to eliminate the current column in row i
				for (int j{ k }; j < 12; ++j) {
					aug[i][j] -= factor * aug[k][j];					// subtract the scaled pivot row from the current row to eliminate the current column
				}
			}
		}

		// Back substitution -- normalize diagonal and eliminate above
		for (int k{ 5 }; k >= 0; --k)
		{
			double pivot{ aug[k][k] };
			for (int j{ 0 }; j < 12; ++j)								// divide entire row by pivot to make leading coefficient 1
				aug[k][j] /= pivot;

																		// eliminate column k from all rows above
			for (int i{ 0 }; i < k; ++i) {
				double factor{ aug[i][k] };
				for (int j{ 0 }; j < 12; ++j)
					aug[i][j] -= factor * aug[k][j];
			}
		}

		return Matrix6x6(
			Vector6(aug[0][6], aug[0][7], aug[0][8], aug[0][9], aug[0][10], aug[0][11]),
			Vector6(aug[1][6], aug[1][7], aug[1][8], aug[1][9], aug[1][10], aug[1][11]),
			Vector6(aug[2][6], aug[2][7], aug[2][8], aug[2][9], aug[2][10], aug[2][11]),
			Vector6(aug[3][6], aug[3][7], aug[3][8], aug[3][9], aug[3][10], aug[3][11]),
			Vector6(aug[4][6], aug[4][7], aug[4][8], aug[4][9], aug[4][10], aug[4][11]),
			Vector6(aug[5][6], aug[5][7], aug[5][8], aug[5][9], aug[5][10], aug[5][11])
		);
	}
};