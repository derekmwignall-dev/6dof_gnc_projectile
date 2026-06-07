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

//// Namespace that contains functions to generate rotation matrices for rotations about the x, y, and z axes by a specified angle in degrees
//namespace Rotation {
//	constexpr Matrix3x3 rotX(double angle) {
//		double rad = angle * Constants::pi / 180.0; // convert to radians
//		double c = std::cos(rad);
//		double s = std::sin(rad);
//		return Matrix3x3(
//			Vector3(1, 0, 0),
//			Vector3(0, c, -s),
//			Vector3(0, s, c)
//		);
//	}
//	constexpr Matrix3x3 rotY(double angle) {
//		double rad = angle * Constants::pi / 180.0; // convert to radians
//		double c = std::cos(rad);
//		double s = std::sin(rad);
//		return Matrix3x3(
//			Vector3(c, 0, s),
//			Vector3(0, 1, 0),
//			Vector3(-s, 0, c)
//		);
//	}
//	constexpr Matrix3x3 rotZ(double angle) {
//		double rad = angle * Constants::pi / 180.0; // convert to radians
//		double c = std::cos(rad);
//		double s = std::sin(rad);
//		return Matrix3x3(
//			Vector3(c, -s, 0),
//			Vector3(s, c, 0),
//			Vector3(0, 0, 1)
//		);
//	}
//}

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