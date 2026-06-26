#pragma once
#include "linalg.h"

// ============================================================
// Matrix3x9 — 3 rows, 9 columns
// ============================================================
class Matrix3x9 {
private:
	Vector9 m_row0{};
	Vector9 m_row1{};
	Vector9 m_row2{};
public:
	constexpr Matrix3x9(const Vector9& row0, const Vector9& row1, const Vector9& row2)
		: m_row0{ row0 }, m_row1{ row1 }, m_row2{ row2 } {
	}
	constexpr Matrix3x9() : m_row0{}, m_row1{}, m_row2{} {} // default to zero matrix

	Vector9 getRow(int i) const
	{
		switch (i) {
		case 0: return m_row0;
		case 1: return m_row1;
		case 2: return m_row2;
		default: throw std::out_of_range("Index out of range for Matrix3x9");
		}
	}

	// (3x9) * M9x9 -> Matrix3x9
	Matrix3x9 operator*(const Matrix9x9& M) const
	{
		// result[i][j] = row_i(this) . col_j(M)
		// col_j(M) = the j-th element of each of M's 9 rows
		auto dotRowCol = [&](const Vector9& row, int col) -> double {
			double sum = 0.0;
			for (int k = 0; k < 9; ++k)
				sum += row.get(k) * M.getRow(k).get(col);
			return sum;
			};

		auto makeRow9 = [&](const Vector9& row) -> Vector9 {
			return Vector9(
				dotRowCol(row, 0), dotRowCol(row, 1), dotRowCol(row, 2),
				dotRowCol(row, 3), dotRowCol(row, 4), dotRowCol(row, 5),
				dotRowCol(row, 6), dotRowCol(row, 7), dotRowCol(row, 8)
			);
			};

		return Matrix3x9(
			makeRow9(getRow(0)),
			makeRow9(getRow(1)),
			makeRow9(getRow(2))
		);
	}

	// H * v ->  Vector3
	Vector3 operator*(const Vector9& v) const
	{
		return Vector3(
			m_row0.dotP(v),
			m_row1.dotP(v),
			m_row2.dotP(v)
		);
	}

	// H * M  ->  Matrix3x3
	Matrix3x3 operator*(const class Matrix9x3& M) const;  // defined after Matrix9x3

	// Hᵀ  ->  Matrix9x3
	Matrix9x3 trans() const;  // defined after Matrix9x3

	constexpr Matrix3x9 operator*(double scalar) const
	{
		return Matrix3x9(m_row0 * scalar, m_row1 * scalar, m_row2 * scalar);
	}

	constexpr Matrix3x9 operator+(const Matrix3x9& u) const
	{
		return Matrix3x9(m_row0 + u.m_row0, m_row1 + u.m_row1, m_row2 + u.m_row2);
	}

	constexpr Matrix3x9 operator-(const Matrix3x9& u) const
	{
		return Matrix3x9(m_row0 - u.m_row0, m_row1 - u.m_row1, m_row2 - u.m_row2);
	}

};


// ============================================================
// Matrix9x3 — 9 rows, 3 columns
// ============================================================
class Matrix9x3 {
private:
	Vector3 m_row0{};
	Vector3 m_row1{};
	Vector3 m_row2{};
	Vector3 m_row3{};
	Vector3 m_row4{};
	Vector3 m_row5{};
	Vector3 m_row6{};
	Vector3 m_row7{};
	Vector3 m_row8{};
public:
	constexpr Matrix9x3(const Vector3& row0, const Vector3& row1, const Vector3& row2,
		const Vector3& row3, const Vector3& row4, const Vector3& row5,
		const Vector3& row6, const Vector3& row7, const Vector3& row8)
		: m_row0{ row0 }, m_row1{ row1 }, m_row2{ row2 },
		m_row3{ row3 }, m_row4{ row4 }, m_row5{ row5 },
		m_row6{ row6 }, m_row7{ row7 }, m_row8{ row8 } {
	}
	constexpr Matrix9x3()
		: m_row0{}, m_row1{}, m_row2{}, m_row3{}, m_row4{},
		m_row5{}, m_row6{}, m_row7{}, m_row8{} {
	} // default to zero matrix

	Vector3 getRow(int i) const
	{
		switch (i) {
		case 0: return m_row0;
		case 1: return m_row1;
		case 2: return m_row2;
		case 3: return m_row3;
		case 4: return m_row4;
		case 5: return m_row5;
		case 6: return m_row6;
		case 7: return m_row7;
		case 8: return m_row8;
		default: throw std::out_of_range("Index out of range for Matrix9x3");
		}
	}

	// H * v  ->  Vector9
	Vector9 operator*(const Vector3& v) const
	{
		return Vector9(
			m_row0.dotP(v),
			m_row1.dotP(v),
			m_row2.dotP(v),
			m_row3.dotP(v),
			m_row4.dotP(v),
			m_row5.dotP(v),
			m_row6.dotP(v),
			m_row7.dotP(v),
			m_row8.dotP(v)
		);
	}

	// H * M ->  Matrix9x9
	Matrix9x9 operator*(const Matrix3x9& H) const
	{
		Matrix3x9 HT_rows{ H }; // we'll dot each of our rows against each column of H
		// column j of H = row j of Hᵀ
		// (K*H)[i][j] = row_i(K) . col_j(H) = row_i(K) . row_j(Hᵀ)
		// Build result row by row
		auto dotRowCol = [&](const Vector3& kRow, int hCol) -> double {
			// column hCol of H = the hCol-th element of each H row
			return kRow.getX() * H.getRow(0).get(hCol)
				+ kRow.getY() * H.getRow(1).get(hCol)
				+ kRow.getZ() * H.getRow(2).get(hCol);
			};

		auto makeRow9 = [&](const Vector3& kRow) -> Vector9 {
			return Vector9(
				dotRowCol(kRow, 0), dotRowCol(kRow, 1), dotRowCol(kRow, 2),
				dotRowCol(kRow, 3), dotRowCol(kRow, 4), dotRowCol(kRow, 5),
				dotRowCol(kRow, 6), dotRowCol(kRow, 7), dotRowCol(kRow, 8)
			);
			};

		return Matrix9x9(
			makeRow9(m_row0), makeRow9(m_row1), makeRow9(m_row2),
			makeRow9(m_row3), makeRow9(m_row4), makeRow9(m_row5),
			makeRow9(m_row6), makeRow9(m_row7), makeRow9(m_row8)
		);
	}

	// Hᵀ  ->  Matrix3x9
	Matrix3x9 trans() const
	{
		return Matrix3x9(
			Vector9(m_row0.getX(), m_row1.getX(), m_row2.getX(), m_row3.getX(), m_row4.getX(), m_row5.getX(), m_row6.getX(), m_row7.getX(), m_row8.getX()),
			Vector9(m_row0.getY(), m_row1.getY(), m_row2.getY(), m_row3.getY(), m_row4.getY(), m_row5.getY(), m_row6.getY(), m_row7.getY(), m_row8.getY()),
			Vector9(m_row0.getZ(), m_row1.getZ(), m_row2.getZ(), m_row3.getZ(), m_row4.getZ(), m_row5.getZ(), m_row6.getZ(), m_row7.getZ(), m_row8.getZ())
		);
	}

	constexpr Matrix9x3 operator*(double scalar) const
	{
		return Matrix9x3(
			m_row0 * scalar, m_row1 * scalar, m_row2 * scalar,
			m_row3 * scalar, m_row4 * scalar, m_row5 * scalar,
			m_row6 * scalar, m_row7 * scalar, m_row8 * scalar
		);
	}

	constexpr Matrix9x3 operator+(const Matrix9x3& u) const
	{
		return Matrix9x3(
			m_row0 + u.m_row0, m_row1 + u.m_row1, m_row2 + u.m_row2,
			m_row3 + u.m_row3, m_row4 + u.m_row4, m_row5 + u.m_row5,
			m_row6 + u.m_row6, m_row7 + u.m_row7, m_row8 + u.m_row8
		);
	}
};


// ============================================================
// Out-of-line definitions requiring both classes to be complete
// ============================================================

// Matrix3x9 * Matrix9x3  ->  Matrix3x3
inline Matrix3x3 Matrix3x9::operator*(const Matrix9x3& M) const
{

	// cleaner: extract each column of M then dot with H rows
	auto col0 = [&](int row) { return M.getRow(row).getX(); };
	auto col1 = [&](int row) { return M.getRow(row).getY(); };
	auto col2 = [&](int row) { return M.getRow(row).getZ(); };

	auto rowDotCol = [&](const Vector9& hRow, auto colFn) -> double {
		double sum{ 0.0 };
		for (int i{ 0 }; i < 9; ++i) sum += hRow.get(i) * colFn(i);
		return sum;
		};

	return Matrix3x3(
		Vector3(rowDotCol(m_row0, col0), rowDotCol(m_row0, col1), rowDotCol(m_row0, col2)),
		Vector3(rowDotCol(m_row1, col0), rowDotCol(m_row1, col1), rowDotCol(m_row1, col2)),
		Vector3(rowDotCol(m_row2, col0), rowDotCol(m_row2, col1), rowDotCol(m_row2, col2))
	);
}

// Matrix3x9::trans()  ->  Matrix9x3
inline Matrix9x3 Matrix3x9::trans() const
{
	// row i of Hᵀ = column i of H = element i of each H row
	return Matrix9x3(
		Vector3(m_row0.get(0), m_row1.get(0), m_row2.get(0)),
		Vector3(m_row0.get(1), m_row1.get(1), m_row2.get(1)),
		Vector3(m_row0.get(2), m_row1.get(2), m_row2.get(2)),
		Vector3(m_row0.get(3), m_row1.get(3), m_row2.get(3)),
		Vector3(m_row0.get(4), m_row1.get(4), m_row2.get(4)),
		Vector3(m_row0.get(5), m_row1.get(5), m_row2.get(5)),
		Vector3(m_row0.get(6), m_row1.get(6), m_row2.get(6)),
		Vector3(m_row0.get(7), m_row1.get(7), m_row2.get(7)),
		Vector3(m_row0.get(8), m_row1.get(8), m_row2.get(8))
	);
}

// Matrix9x9 * Matrix9x3  ->  Matrix9x3
inline Matrix9x3 operator*(const Matrix9x9& A, const Matrix9x3& B)
{
	auto col0 = [&](int row) { return B.getRow(row).getX(); };
	auto col1 = [&](int row) { return B.getRow(row).getY(); };
	auto col2 = [&](int row) { return B.getRow(row).getZ(); };

	auto rowDotCol = [&](const Vector9& aRow, auto colFn) -> double {
		double sum{ 0.0 };
		for (int i{ 0 }; i < 9; ++i) sum += aRow.get(i) * colFn(i);
		return sum;
		};

	auto makeRow = [&](int r) -> Vector3 {
		Vector9 aRow{ A.getRow(r) };
		return Vector3(rowDotCol(aRow, col0), rowDotCol(aRow, col1), rowDotCol(aRow, col2));
		};

	return Matrix9x3(
		makeRow(0), makeRow(1), makeRow(2),
		makeRow(3), makeRow(4), makeRow(5),
		makeRow(6), makeRow(7), makeRow(8)
	);
}

// Matrix9x3 * Matrix3x3  ->  Matrix9x3
inline Matrix9x3 operator*(const Matrix9x3& A, const Matrix3x3& B)
{
	Matrix3x3 BT{ B.trans() };

	auto makeRow = [&](const Vector3& aRow) -> Vector3 {
		return Vector3(
			aRow.dotP(BT.getRow(0)),
			aRow.dotP(BT.getRow(1)),
			aRow.dotP(BT.getRow(2))
		);
		};

	return Matrix9x3(
		makeRow(A.getRow(0)), makeRow(A.getRow(1)), makeRow(A.getRow(2)),
		makeRow(A.getRow(3)), makeRow(A.getRow(4)), makeRow(A.getRow(5)),
		makeRow(A.getRow(6)), makeRow(A.getRow(7)), makeRow(A.getRow(8))
	);
}