#pragma once

#undef min
#undef max

#include <cmath>

struct vec2; struct vec3; struct vec4;
struct mat2; struct mat3; struct mat4;
struct quat;

inline float pi_constant() { return 3.1415927410125732421875f; }
inline float degrees(float rad) { return rad * 180.0f / pi_constant(); }
inline float radians(float deg) { return deg * pi_constant() / 180.0f; }
inline float clamp(float v, float min, float max) { return v < min ? min : v > max ? max : v; }
inline float mix(float a, float b, float t) { return a + (b - a) * t; }

#pragma mark - vec2

struct vec2
{
	float x, y;

	vec2() {}
	vec2(vec2 const& v) : x(v.x), y(v.y) {}
	explicit vec2(float s) : x(s), y(s) {}
	explicit vec2(float x, float y) : x(x), y(y) {}
	explicit vec2(vec3 const& v);

	float& operator [] (std::size_t i) { return (&x)[i]; }
	float const& operator [] (std::size_t i) const { return (&x)[i]; }

	vec2& operator += (float s) { x += s; y += s; return *this; }
	vec2& operator -= (float s) { x -= s; y -= s; return *this; }
	vec2& operator *= (float s) { x *= s; y *= s; return *this; }
	vec2& operator /= (float s) { float const i = 1 / s; x *= i; y *= i; return *this; }
	vec2& operator += (vec2 const& v) { x += v.x; y += v.y; return *this; }
	vec2& operator -= (vec2 const& v) { x -= v.x; y -= v.y; return *this; }
	vec2& operator *= (vec2 const& v) { x *= v.x; y *= v.y; return *this; }
	vec2& operator /= (vec2 const& v) { x /= v.x; y /= v.y; return *this; }
};

inline vec2  operator + (vec2 const& v, float s) { return vec2(v.x + s, v.y + s); }
inline vec2  operator - (vec2 const& v, float s) { return vec2(v.x - s, v.y - s); }
inline vec2  operator * (vec2 const& v, float s) { return vec2(v.x * s, v.y * s); }
inline vec2  operator / (vec2 const& v, float s) { float const i = 1 / s; return vec2(v.x * i, v.y * i); }
inline vec2  operator + (float s, vec2 const& v) { return vec2(s + v.x, s + v.y); }
inline vec2  operator - (float s, vec2 const& v) { return vec2(s - v.x, s - v.y); }
inline vec2  operator * (float s, vec2 const& v) { return vec2(s * v.x, s * v.y); }
inline vec2  operator / (float s, vec2 const& v) { return vec2(s / v.x, s / v.y); }
inline vec2  operator + (vec2 const& a, vec2 const& b) { return vec2(a.x + b.x, a.y + b.y); }
inline vec2  operator - (vec2 const& a, vec2 const& b) { return vec2(a.x - b.x, a.y - b.y); }
inline vec2  operator * (vec2 const& a, vec2 const& b) { return vec2(a.x * b.x, a.y * b.y); }
inline vec2  operator / (vec2 const& a, vec2 const& b) { return vec2(a.x / b.x, a.y / b.y); }
inline vec2  inverse(vec2 const& v) { return vec2(-v.x, -v.y); }
inline vec2  abs(vec2 const& v) { return vec2(std::abs(v.x), std::abs(v.y)); }
inline vec2  ceil(vec2 const& v) { return vec2(std::ceil(v.x), std::ceil(v.y)); }
inline vec2  floor(vec2 const& v) { return vec2(std::floor(v.x), std::floor(v.y)); }
inline vec2  min(vec2 const& a, vec2 const& b) { return vec2(std::fmin(a.x, b.x), std::fmin(a.y, b.y)); }
inline vec2  max(vec2 const& a, vec2 const& b) { return vec2(std::fmax(a.x, b.x), std::fmax(a.y, b.y)); }
inline vec2  clamp(vec2 const& v, vec2 const& min, vec2 const& max) { return vec2(clamp(v.x, min.x, max.x), clamp(v.y, min.y, max.y)); }
inline float dot(vec2 const& a, vec2 const& b) { return a.x * b.x + a.y * b.y; }
inline float length(vec2 const& v) { return std::sqrt(dot(v, v)); }
inline float distance(vec2 const& a, vec2 const& b) { return length(a - b); }
inline vec2  normalize(vec2 const& v) { return v / length(v); }
inline vec2  mix(vec2 const& a, vec2 const& b, float t) { return vec2(mix(a.x, b.x, t), mix(a.y, b.y, t)); }
inline vec2  reflect(vec2 const& i, vec2 const& n) { return i - n * dot(n, i) * 2; }

inline vec2 refract(vec2 const& i, vec2 const& n, float eta)
{
	float const dni = dot(n, i);
	float const k = 1 - eta * eta * (1 - dni * dni);
	return k < 0 ? vec2(0) : (i * eta - n * (eta * dni + std::sqrt(k)));
}

#pragma mark - vec3

struct vec3
{
	float x, y, z;

	vec3() {}
	vec3(vec3 const& v) : x(v.x), y(v.y), z(v.z) {}
	explicit vec3(float s) : x(s), y(s), z(s) {}
	explicit vec3(float x, float y, float z) : x(x), y(y), z(z) {}
	explicit vec3(vec2 const& v, float z) : x(v.x), y(v.y), z(z) {}
	explicit vec3(vec4 const&);

	float& operator [] (std::size_t i) { return (&x)[i]; }
	float const& operator [] (std::size_t i) const { return (&x)[i]; }

	vec3& operator += (float s) { x += s; y += s; z += s; return *this; }
	vec3& operator -= (float s) { x -= s; y -= s; z -= s; return *this; }
	vec3& operator *= (float s) { x *= s; y *= s; z *= s; return *this; }
	vec3& operator /= (float s) { float const i = 1 / s; x *= i; y *= i; z *= i; return *this; }
	vec3& operator += (vec3 const& v) { x += v.x; y += v.y; z += v.z; return *this; }
	vec3& operator -= (vec3 const& v) { x -= v.x; y -= v.y; z -= v.z; return *this; }
	vec3& operator *= (vec3 const& v) { x *= v.x; y *= v.y; z *= v.z; return *this; }
	vec3& operator /= (vec3 const& v) { x /= v.x; y /= v.y; z /= v.z; return *this; }
};

inline vec3  operator + (vec3 const& v, float s) { return vec3(v.x + s, v.y + s, v.z + s); }
inline vec3  operator - (vec3 const& v, float s) { return vec3(v.x - s, v.y - s, v.z - s); }
inline vec3  operator * (vec3 const& v, float s) { return vec3(v.x * s, v.y * s, v.z * s); }
inline vec3  operator / (vec3 const& v, float s) { float const i = 1 / s; return vec3(v.x * i, v.y * i, v.z * i); }
inline vec3  operator + (float s, vec3 const& v) { return vec3(s + v.x, s + v.y, s + v.z); }
inline vec3  operator - (float s, vec3 const& v) { return vec3(s - v.x, s - v.y, s - v.z); }
inline vec3  operator * (float s, vec3 const& v) { return vec3(s * v.x, s * v.y, s * v.z); }
inline vec3  operator / (float s, vec3 const& v) { return vec3(s / v.x, s / v.y, s / v.z); }
inline vec3  operator + (vec3 const& a, vec3 const& b) { return vec3(a.x + b.x, a.y + b.y, a.z + b.z); }
inline vec3  operator - (vec3 const& a, vec3 const& b) { return vec3(a.x - b.x, a.y - b.y, a.z - b.z); }
inline vec3  operator * (vec3 const& a, vec3 const& b) { return vec3(a.x * b.x, a.y * b.y, a.z * b.z); }
inline vec3  operator / (vec3 const& a, vec3 const& b) { return vec3(a.x / b.x, a.y / b.y, a.z / b.z); }
inline vec3  inverse(vec3 const& v) { return vec3(-v.x, -v.y, -v.z); }
inline vec3  min(vec3 const& a, vec3 const& b) { return vec3(std::fmin(a.x, b.x), std::fmin(a.y, b.y), std::fmin(a.z, b.z)); }
inline vec3  max(vec3 const& a, vec3 const& b) { return vec3(std::fmax(a.x, b.x), std::fmax(a.y, b.y), std::fmax(a.z, b.z)); }
inline vec3  clamp(vec3 const& v, vec3 const& min, vec3 const& max) { return vec3(clamp(v.x, min.x, max.x), clamp(v.y, min.y, max.y), clamp(v.z, min.z, max.z)); }
inline vec3  abs(vec3 const& v) { return vec3(std::abs(v.x), std::abs(v.y), std::abs(v.z)); }
inline vec3  floor(vec3 const& v) { return vec3(std::floor(v.x), std::floor(v.y), std::floor(v.z)); }
inline vec3  ceil(vec3 const& v) { return vec3(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z)); }
inline float dot(vec3 const& a, vec3 const& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline float length(vec3 const& v) { return std::sqrt(dot(v, v)); }
inline float distance(vec3 const& a, vec3 const& b) { return length(a - b); }
inline vec3  normalize(vec3 const& v) { return v / length(v); }
inline vec3  cross(vec3 const& a, vec3 const& b) { return vec3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x); }
inline vec3  mix(vec3 const& a, vec3 const& b, float t) { return vec3(mix(a.x, b.x, t), mix(a.y, b.y, t), mix(a.z, b.z, t)); }
inline vec3  reflect(vec3 const& i, vec3 const& n) { return i - n * dot(n, i) * 2; }

inline vec3 refract(vec3 const& i, vec3 const& n, float eta)
{
	float const dni = dot(n, i);
	float const k = 1 - eta * eta * (1 - dni * dni);
	return k < 0 ? vec3(0) : (i * eta - n * (eta * dni + std::sqrt(k)));
}

#pragma mark - vec4

struct vec4
{
	float x, y, z, w;

	vec4() {}
	vec4(vec4 const& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}
	explicit vec4(float s) : x(s), y(s), z(s), w(s) {}
	explicit vec4(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
	explicit vec4(vec2 const& v, float z, float w) : x(v.x), y(v.y), z(z), w(w) {}
	explicit vec4(vec3 const& v, float w) : x(v.x), y(v.y), z(v.z), w(w) {}

	float& operator [] (std::size_t i) { return (&x)[i]; }
	float const& operator [] (std::size_t i) const { return (&x)[i]; }

	vec4& operator += (float s) { x += s; y += s; z += s; w += s; return *this; }
	vec4& operator -= (float s) { x -= s; y -= s; z -= s; w -= s; return *this; }
	vec4& operator *= (float s) { x *= s; y *= s; z *= s; w *= s; return *this; }
	vec4& operator /= (float s) { float const i = 1 / s; x *= i; y *= i; z *= i; w *= i; return *this; }
	vec4& operator += (vec4 const& v) { x += v.x; y += v.y; z += v.z; w += v.w; return *this; }
	vec4& operator -= (vec4 const& v) { x -= v.x; y -= v.y; z -= v.z; w -= v.w; return *this; }
	vec4& operator *= (vec4 const& v) { x *= v.x; y *= v.y; z *= v.z; w *= v.w; return *this; }
	vec4& operator /= (vec4 const& v) { x /= v.x; y /= v.y; z /= v.z; w /= v.w; return *this; }
};

inline vec4  operator + (vec4 const& v, float s) { return vec4(v.x + s, v.y + s, v.z + s, v.w + s); }
inline vec4  operator - (vec4 const& v, float s) { return vec4(v.x - s, v.y - s, v.z - s, v.w - s); }
inline vec4  operator * (vec4 const& v, float s) { return vec4(v.x * s, v.y * s, v.z * s, v.w * s); }
inline vec4  operator / (vec4 const& v, float s) { float const i = 1 / s; return vec4(v.x * i, v.y * i, v.z * i, v.w * i); }
inline vec4  operator + (float s, vec4 const& v) { return vec4(s + v.x, s + v.y, s + v.z, s + v.w); }
inline vec4  operator - (float s, vec4 const& v) { return vec4(s - v.x, s - v.y, s - v.z, s - v.w); }
inline vec4  operator * (float s, vec4 const& v) { return vec4(s * v.x, s * v.y, s * v.z, s * v.w); }
inline vec4  operator / (float s, vec4 const& v) { return vec4(s / v.x, s / v.y, s / v.z, s / v.w); }
inline vec4  operator + (vec4 const& a, vec4 const& b) { return vec4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w); }
inline vec4  operator - (vec4 const& a, vec4 const& b) { return vec4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w); }
inline vec4  operator * (vec4 const& a, vec4 const& b) { return vec4(a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w); }
inline vec4  operator / (vec4 const& a, vec4 const& b) { return vec4(a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w); }
inline vec4  inverse(vec4 const& v) { return vec4(-v.x, -v.y, -v.z, -v.w); }
inline vec4  min(vec4 const& a, vec4 const& b) { return vec4(std::fmin(a.x, b.x), std::fmin(a.y, b.y), std::fmin(a.z, b.z), std::fmin(a.w, b.w)); }
inline vec4  max(vec4 const& a, vec4 const& b) { return vec4(std::fmax(a.x, b.x), std::fmax(a.y, b.y), std::fmax(a.z, b.z), std::fmax(a.w, b.w)); }
inline vec4  clamp(vec4 const& v, vec4 const& min, vec4 const& max) { return vec4(clamp(v.x, min.x, max.x), clamp(v.y, min.y, max.y), clamp(v.z, min.z, max.z), clamp(v.w, min.w, max.w)); }
inline vec4  abs(vec4 const& v) { return vec4(std::abs(v.x), std::abs(v.y), std::abs(v.z), std::abs(v.w)); }
inline vec4  floor(vec4 const& v) { return vec4(std::floor(v.x), std::floor(v.y), std::floor(v.z), std::floor(v.w)); }
inline vec4  ceil(vec4 const& v) { return vec4(std::ceil(v.x), std::ceil(v.y), std::ceil(v.z), std::ceil(v.w)); }
inline float dot(vec4 const& a, vec4 const& b) { return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w; }
inline float length(vec4 const& v) { return std::sqrt(dot(v, v)); }
inline float distance(vec4 const& a, vec4 const& b) { return length(a - b); }
inline vec4  normalize(vec4 const& v) { return v / length(v); }
inline vec4  mix(vec4 const& a, vec4 const& b, float t) { return vec4(mix(a.x, b.x, t), mix(a.y, b.y, t), mix(a.z, b.z, t), mix(a.w, b.w, t)); }
inline vec4  reflect(vec4 const& i, vec4 const& n) { return i - n * dot(n, i) * 2; }

inline vec4 refract(vec4 const& i, vec4 const& n, float eta)
{
	float const dni = dot(n, i);
	float const k = 1 - eta * eta * (1 - dni * dni);
	return k < 0 ? vec4(0) : (i * eta - n * (eta * dni + std::sqrt(k)));
}

#pragma mark - mat2

struct mat2
{
	vec2 x, y;

	mat2() {}
	explicit mat2(vec2 const& x, vec2 const& y) : x(x), y(y) {}

	vec2& operator [] (std::size_t i) { return (&x)[i]; }
	vec2 const& operator [] (std::size_t i) const { return (&x)[i]; }

	float determinant() const { return x.x * y.y - x.y * y.x; }

	static mat2 const& identity()
	{
		static mat2 mat(vec2(1, 0), vec2(0, 1));
		return mat;
	}

	static mat2 from_angle(float angle)
	{
		float sa = std::sin(angle);
		float ca = std::cos(angle);
		return mat2(vec2(ca, sa), vec2(-sa, ca));
	}
};

inline mat2 transpose(mat2 const& m)
{
	return mat2(vec2(m.x.x, m.y.x), vec2(m.x.y, m.y.y));
}

inline mat2 inverse(mat2 const& m)
{
	float const d = 1 / m.determinant();
	if (d != 0)
	{
		float const id = 1 / d;
		return mat2(
			vec2(m.y.y, -m.x.y) * id,
			vec2(-m.y.x, m.x.x) * id);
	}
	else
	{
		return mat2::identity();
	}
}

inline vec2 operator * (mat2 const& m, vec2 const& v)
{
	return vec2(
		v.x * m.x.x + v.y * m.y.x,
		v.x * m.x.y + v.y * m.y.y);
}

inline mat2 operator * (mat2 const& a, mat2 const& b)
{
	return mat2(
		vec2(
		a.x.x * b.x.x + a.x.y * b.y.x,
		a.x.x * b.x.y + a.x.y * b.y.y),
		vec2(
		a.y.x * b.x.x + a.y.y * b.y.x,
		a.y.x * b.x.y + a.y.y * b.y.y));
}

#pragma mark - mat3

struct mat3
{
	vec3 x, y, z;

	mat3() {}
	explicit mat3(vec3 const& x, vec3 const& y, vec3 const& z) : x(x), y(y), z(z) {}
	explicit mat3(mat2 const& m, vec3 const& z) : x(m.x, 0), y(m.y, 0), z(z) {}
	explicit mat3(mat4 const&);
	explicit mat3(quat const&);

	vec3& operator [] (std::size_t i) { return (&x)[i]; }
	vec3 const& operator [] (std::size_t i) const { return (&x)[i]; }

	float operator () (int r, int c) const
	{
		return (*this)[c][r];
	}

	float determinant() const
	{
		return
			x.x * (y.y * z.z - y.z * z.y) -
			x.y * (y.x * z.z - y.z * z.x) +
			x.z * (y.x * z.y - y.y * z.x);
	}

	static mat3 const& identity()
	{
		static mat3 mat(vec3(1, 0, 0), vec3(0, 1, 0), vec3(0, 0, 1));
		return mat;
	}

	static mat3 from_axis_angle(vec3 const& axis, float angle)
	{
		float const xy = axis.x * axis.y;
		float const xz = axis.x * axis.z;
		float const yz = axis.y * axis.z;

		float s = std::sin(angle);
		float c = 1 - std::cos(angle);

		return mat3(
			vec3(
			1 + c * (axis.x * axis.x - 1),
			-axis.z * s + c * xy,
			axis.y * s + c * xz),
			vec3(
			axis.z * s + c * xy,
			1 + c * (axis.y * axis.y - 1),
			-axis.x * s + c * yz),
			vec3(
			-axis.y * s + c * xz,
			axis.x * s + c * yz,
			1 + c * (axis.z * axis.z - 1)));
	}

	static mat3 from_euler_angles(float x, float y, float z)
	{
		float sx = std::sin(x); float cx = std::cos(x);
		float sy = std::sin(y); float cy = std::cos(y);
		float sz = std::sin(z); float cz = std::cos(z);
		return mat3(
			vec3(cy * cz, sy * sx - cy * sz * cx, cy * sz * sx + sy * cx),
			vec3(sz, cz * cx, -cz * sx),
			vec3(-sy * cz, sy * sz * cx + cy * sx, cy * cx - sy * sz * sx));
	}
};

inline mat3 transpose(mat3 const& m)
{
	return mat3(
		vec3(m.x.x, m.y.x, m.z.x),
		vec3(m.x.y, m.y.y, m.z.y),
		vec3(m.x.z, m.y.z, m.z.z));
}

inline mat3 inverse(mat3 const& m)
{
	float const d = 1 / m.determinant();
	if (d != 0)
	{
		float const id = 1 / d;
		return mat3(
			vec3(
			id * (m.y.y * m.z.z - m.y.z * m.z.y),
			-id * (m.x.y * m.z.z - m.x.z * m.z.y),
			id * (m.x.y * m.y.z - m.x.z * m.y.y)),

			vec3(
			-id * (m.y.x * m.z.z - m.y.z * m.z.x),
			id * (m.x.x * m.z.z - m.x.z * m.z.x),
			-id * (m.x.x * m.y.z - m.x.z * m.y.x)),

			vec3(
			id * (m.y.x * m.z.y - m.y.y * m.z.x),
			-id * (m.x.x * m.z.y - m.x.y * m.z.x),
			id * (m.x.x * m.y.y - m.x.y * m.y.x)));
	}
	else
	{
		return mat3::identity();
	}
}

inline vec3 operator * (mat3 const& m, vec3 const& v)
{
	return vec3(
		v.x * m.x.x + v.y * m.y.x + v.z * m.z.x,
		v.x * m.x.y + v.y * m.y.y + v.z * m.z.y,
		v.x * m.x.z + v.y * m.y.z + v.z * m.z.z);
}

inline mat3 operator * (mat3 const& a, mat3 const& b)
{
	return mat3(
		vec3(
		a.x.x * b.x.x + a.x.y * b.y.x + a.x.z * b.z.x,
		a.x.x * b.x.y + a.x.y * b.y.y + a.x.z * b.z.y,
		a.x.x * b.x.z + a.x.y * b.y.z + a.x.z * b.z.z),
		vec3(
		a.y.x * b.x.x + a.y.y * b.y.x + a.y.z * b.z.x,
		a.y.x * b.x.y + a.y.y * b.y.y + a.y.z * b.z.y,
		a.y.x * b.x.z + a.y.y * b.y.z + a.y.z * b.z.z),
		vec3(
		a.z.x * b.x.x + a.z.y * b.y.x + a.z.z * b.z.x,
		a.z.x * b.x.y + a.z.y * b.y.y + a.z.z * b.z.y,
		a.z.x * b.x.z + a.z.y * b.y.z + a.z.z * b.z.z));
}

#pragma mark - mat4

mat4 operator * (mat4 const& a, mat4 const& b);

struct mat4
{
	vec4 x, y, z, w;

	mat4() {}
	explicit mat4(vec4 const& x, vec4 const& y, vec4 const& z, vec4 const& w) : x(x), y(y), z(z), w(w) {}
	explicit mat4(mat3 const& m, vec4 const& w) : x(m.x, 0), y(m.y, 0), z(m.z, 0), w(w) {}

	vec4& operator [] (std::size_t i) { return (&x)[i]; }
	vec4 const& operator [] (std::size_t i) const { return (&x)[i]; }

	static mat4 const& identity()
	{
		static mat4 mat(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(0, 0, 0, 1));
		return mat;
	}

	static mat4 frustum(float left, float right, float bottom, float top, float znear, float zfar)
	{
		float const a = 2 * znear;
		float const b = right - left;
		float const c = top - bottom;
		float const d = zfar - znear;
		return mat4(
			vec4(a / b, 0, 0, 0),
			vec4(0, a / c, 0, 0),
			vec4((right + left) / b, (top + bottom) / c, (-zfar - znear) / d, -1),
			vec4(0, 0, (-a * zfar) / d, 0));
	}

	static mat4 ortho(float left, float right, float bottom, float top, float znear, float zfar)
	{
		float const rl = right - left;
		float const tb = top - bottom;
		float const fn = zfar - znear;

		return mat4(
			vec4(2 / rl, 0, 0, 0),
			vec4(0, 2 / tb, 0, 0),
			vec4(0, 0, -2 / fn, 0),
			vec4(-(right + left) / rl, -(top + bottom) / tb, -(zfar + znear) / fn, 1));
	}

	static mat4 perspective(float width, float height, float fov_in_radians, float znear, float zfar)
	{
		float const ymax = znear * std::tan(fov_in_radians / 2);
		float const xmax = ymax * (width / height);
		return frustum(-xmax, xmax, -ymax, ymax, znear, zfar);
	}

	static mat4 look_at(vec3 const& eye, vec3 const& target, vec3 const& up)
	{
		// TODO: clean up
		vec3 const z = normalize(eye - target);
		vec3 const x = normalize(cross(up, z));
		vec3 const y = normalize(cross(z, x));

		mat4 r = mat4(vec4(x.x, y.x, z.x, 0), vec4(x.y, y.y, z.y, 0), vec4(x.z, y.z, z.z, 0), vec4(0, 0, 0, 1));
		mat4 t = mat4(vec4(1, 0, 0, 0), vec4(0, 1, 0, 0), vec4(0, 0, 1, 0), vec4(eye * -1, 1));

		return t * r;
	}
};

inline mat4 transpose(mat4 const& m)
{
	return mat4(
		vec4(m.x.x, m.y.x, m.z.x, m.w.x),
		vec4(m.x.y, m.y.y, m.z.y, m.w.y),
		vec4(m.x.z, m.y.z, m.z.z, m.w.z),
		vec4(m.x.w, m.y.w, m.z.w, m.w.w));
}

inline mat4 inverse(mat4 const& m)
{
	float k[24] =
	{
		m.z.z * m.w.w, m.w.z * m.z.w, m.y.z * m.w.w, m.w.z * m.y.w,
		m.y.z * m.z.w, m.z.z * m.y.w, m.x.z * m.w.w, m.w.z * m.x.w,
		m.x.z * m.z.w, m.z.z * m.x.w, m.x.z * m.y.w, m.y.z * m.x.w,
		m.z.x * m.w.y, m.w.x * m.z.y, m.y.x * m.w.y, m.w.x * m.y.y,
		m.y.x * m.z.y, m.z.x * m.y.y, m.x.x * m.w.y, m.w.x * m.x.y,
		m.x.x * m.z.y, m.z.x * m.x.y, m.x.x * m.y.y, m.y.x * m.x.y
	};

	mat4 r;

	r.x.x = (k[0] * m.y.y + k[3] * m.z.y + k[4] * m.w.y) - (k[1] * m.y.y + k[2] * m.z.y + k[5] * m.w.y);
	r.x.y = (k[1] * m.x.y + k[6] * m.z.y + k[9] * m.w.y) - (k[0] * m.x.y + k[7] * m.z.y + k[8] * m.w.y);
	r.x.z = (k[2] * m.x.y + k[7] * m.y.y + k[10] * m.w.y) - (k[3] * m.x.y + k[6] * m.y.y + k[11] * m.w.y);
	r.x.w = (k[5] * m.x.y + k[8] * m.y.y + k[11] * m.z.y) - (k[4] * m.x.y + k[9] * m.y.y + k[10] * m.z.y);

	r.y.x = (k[1] * m.y.x + k[2] * m.z.x + k[5] * m.w.x) - (k[0] * m.y.x + k[3] * m.z.x + k[4] * m.w.x);
	r.y.y = (k[0] * m.x.x + k[7] * m.z.x + k[8] * m.w.x) - (k[1] * m.x.x + k[6] * m.z.x + k[9] * m.w.x);
	r.y.z = (k[3] * m.x.x + k[6] * m.y.x + k[11] * m.w.x) - (k[2] * m.x.x + k[7] * m.y.x + k[10] * m.w.x);
	r.y.w = (k[4] * m.x.x + k[9] * m.y.x + k[10] * m.z.x) - (k[5] * m.x.x + k[8] * m.y.x + k[11] * m.z.x);

	r.z.x = (k[12] * m.y.w + k[15] * m.z.w + k[16] * m.w.w) - (k[13] * m.y.w + k[14] * m.z.w + k[17] * m.w.w);
	r.z.y = (k[13] * m.x.w + k[18] * m.z.w + k[21] * m.w.w) - (k[12] * m.x.w + k[19] * m.z.w + k[20] * m.w.w);
	r.z.z = (k[14] * m.x.w + k[19] * m.y.w + k[22] * m.w.w) - (k[15] * m.x.w + k[18] * m.y.w + k[23] * m.w.w);
	r.z.w = (k[17] * m.x.w + k[20] * m.y.w + k[23] * m.z.w) - (k[16] * m.x.w + k[21] * m.y.w + k[22] * m.z.w);

	r.w.x = (k[14] * m.z.z + k[17] * m.w.z + k[13] * m.y.z) - (k[16] * m.w.z + k[12] * m.y.z + k[15] * m.z.z);
	r.w.y = (k[20] * m.w.z + k[12] * m.x.z + k[19] * m.z.z) - (k[18] * m.z.z + k[21] * m.w.z + k[13] * m.x.z);
	r.w.z = (k[18] * m.y.z + k[23] * m.w.z + k[15] * m.x.z) - (k[22] * m.w.z + k[14] * m.x.z + k[19] * m.y.z);
	r.w.w = (k[22] * m.z.z + k[16] * m.x.z + k[21] * m.y.z) - (k[20] * m.y.z + k[23] * m.z.z + k[17] * m.x.z);

	float const det = m.x.x * r.x.x + m.y.x * r.x.y + m.z.x * r.x.z + m.w.x * r.x.w;

	if (det != 0)
	{
		float const idet = 1 / det;
		r.x *= idet;
		r.y *= idet;
		r.z *= idet;
		r.w *= idet;
		return r;
	}
	else
	{
		return mat4::identity();
	}
}

inline vec4 operator * (mat4 const& m, vec4 const& v)
{
	return vec4(
		v.x * m.x.x + v.y * m.y.x + v.z * m.z.x + v.w * m.w.x,
		v.x * m.x.y + v.y * m.y.y + v.z * m.z.y + v.w * m.w.y,
		v.x * m.x.z + v.y * m.y.z + v.z * m.z.z + v.w * m.w.z,
		v.x * m.x.w + v.y * m.y.w + v.z * m.z.w + v.w * m.w.w);
}

inline vec3 operator * (mat4 const& m, vec3 const& v)
{
	return vec3(m * vec4(v, 1));
}

inline mat4 operator * (mat4 const& a, mat4 const& b)
{
	return mat4(
		vec4(
		a.x.x * b.x.x + a.x.y * b.y.x + a.x.z * b.z.x + a.x.w * b.w.x,
		a.x.x * b.x.y + a.x.y * b.y.y + a.x.z * b.z.y + a.x.w * b.w.y,
		a.x.x * b.x.z + a.x.y * b.y.z + a.x.z * b.z.z + a.x.w * b.w.z,
		a.x.x * b.x.w + a.x.y * b.y.w + a.x.z * b.z.w + a.x.w * b.w.w),
		vec4(
		a.y.x * b.x.x + a.y.y * b.y.x + a.y.z * b.z.x + a.y.w * b.w.x,
		a.y.x * b.x.y + a.y.y * b.y.y + a.y.z * b.z.y + a.y.w * b.w.y,
		a.y.x * b.x.z + a.y.y * b.y.z + a.y.z * b.z.z + a.y.w * b.w.z,
		a.y.x * b.x.w + a.y.y * b.y.w + a.y.z * b.z.w + a.y.w * b.w.w),
		vec4(
		a.z.x * b.x.x + a.z.y * b.y.x + a.z.z * b.z.x + a.z.w * b.w.x,
		a.z.x * b.x.y + a.z.y * b.y.y + a.z.z * b.z.y + a.z.w * b.w.y,
		a.z.x * b.x.z + a.z.y * b.y.z + a.z.z * b.z.z + a.z.w * b.w.z,
		a.z.x * b.x.w + a.z.y * b.y.w + a.z.z * b.z.w + a.z.w * b.w.w),
		vec4(
		a.w.x * b.x.x + a.w.y * b.y.x + a.w.z * b.z.x + a.w.w * b.w.x,
		a.w.x * b.x.y + a.w.y * b.y.y + a.w.z * b.z.y + a.w.w * b.w.y,
		a.w.x * b.x.z + a.w.y * b.y.z + a.w.z * b.z.z + a.w.w * b.w.z,
		a.w.x * b.x.w + a.w.y * b.y.w + a.w.z * b.z.w + a.w.w * b.w.w));
}

#pragma mark - quat

quat operator * (quat const& a, quat const& b);
quat normalize(quat const& q);

struct quat
{
	float x, y, z, w;

	quat() {}
	explicit quat(float x, float y, float z, float w) : x(x), y(y), z(z), w(w) {}
	explicit quat(vec3 const& v, float s) : x(v.x), y(v.y), z(v.z), w(s) {}

	explicit quat(mat3 const& m) // XXX: needs to be checked
	{
		float tr = m[0][0] + m[1][1] + m[2][2], h;
		if (tr >= 0)
		{
			h = sqrtf(tr + 1);
			w = 0.5f * h;
			h = 0.5f / h;

			x = (m[1][2] - m[2][1]) * h;
			y = (m[2][0] - m[0][2]) * h;
			z = (m[0][1] - m[1][0]) * h;
		}
		else
		{
			unsigned int i = 0;
			if (m[1][1] > m[0][0]) i = 1;
			if (m[2][2] > m[i][i]) i = 2;
			int j = (i + 1) % 3;
			int k = (i + 2) % 3;

			float* q = &x;
			h = sqrtf((m[i][i] - (m[j][j] + m[k][k])) + 1.0);
			q[i] = 0.5f * h;
			h = 0.5f / h;
			q[j] = (m[j][i] + m[i][j]) * h;
			q[k] = (m[i][k] + m[k][i]) * h;
			w = (m[j][k] - m[k][j]) * h;
		}
	}

	quat& operator *= (quat const& q) { return *this = *this * q; }

	vec3 vector() const { return vec3(x, y, z); }
	float norm() const { return std::sqrt(x * x + y * y + z * z + w * w); }
	static quat const& identity()
	{
		static quat q(0, 0, 0, 1);
		return q;
	}

	void to_euler_angles(float& ax, float& ay, float& az)
	{
		float x2 = x * x;
		float y2 = y * y;
		float z2 = z * z;
		float w2 = w * w;

		az = atan2(2.0 * (x * y + z * w), (x2 - y2 - z2 + w2));
		ax = atan2(2.0 * (y * z + x * w), (-x2 - y2 + z2 + w2));
		ay = asin(-2.0 * (x * z - y * w));
	}

	static quat from_euler_angles(float x, float y, float z)
	{
		x *= 0.5; y *= 0.5; z *= 0.5;
		float sx = std::sin(x); float cx = std::cos(x);
		float sy = std::sin(y); float cy = std::cos(y);
		float sz = std::sin(z); float cz = std::cos(z);
		return quat(
			cz * sy * cx + sz * cy * sx,
			cz * cy * sx - sz * sy * cx,
			sz * cy * cx - cz * sy * sx,
			cz * cy * cx + sz * sy * sx);
	}

	static quat from_axis_angle(vec3 const& axis, float angle)
	{
		float len = length(axis);
		if (len == 0)
		{
			return quat::identity();
		}
		else
		{
			angle *= 0.5;
			float s = std::sin(angle);
			float c = std::cos(angle);
			return quat(axis.x * s, axis.y * s, axis.z * s, c);
		}
	}

	static quat from_shortest_arc(vec3 const& a, vec3 const& b)
	{
		vec3 c = cross(a, b);
		quat q(c.x, c.y, c.z, dot(a, b));
		q = normalize(q);
		q.w += 1;
		q = normalize(q);
		return q;
	}
};

inline quat operator * (quat const& a, quat const& b)
{
	return quat(
		a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
		a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z,
		a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.x,
		a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z);
}

inline quat normalize(quat const& q)
{
	float n = q.norm();
	float in = n == 0 ? 1 : 1 / n;
	return quat(q.x * in, q.y * in, q.z * in, q.w * in);
}

inline quat conjugate(quat const& q)
{
	return quat(-q.x, -q.y, -q.z, q.w);
}

inline quat inverse(quat const& q)
{
	quat r = conjugate(q);
	float in = 1 / q.norm();
	return quat(r.x * in, r.y * in, r.z * in, r.w * in);
}

inline vec3 rotate(vec3 const& v, quat const& q)
{
	return ((q * quat(v, 0.0)) * conjugate(q)).vector();
}

inline quat slerp(quat const& a, quat const& b, float t)
{
	float const epsilon = static_cast<float>(1.0e-8);

	float cosine = a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
	float sine = 1 - cosine * cosine;

	float sign;
	if (cosine < 0)
	{
		cosine = -cosine;
		sign = -1;
	}
	else
	{
		sign = 1;
	}

	if (sine >= epsilon * epsilon)
	{
		sine = std::sqrt(sine);

		float const angle = std::atan2(sine, cosine);
		float const i_sin_angle = 1 / sine;

		float lower_weight = std::sin(angle * (1 - t)) * i_sin_angle;
		float upper_weight = std::sin(angle * t) * i_sin_angle * sign;

		return quat(
			a.x * lower_weight + b.x * upper_weight,
			a.y * lower_weight + b.y * upper_weight,
			a.z * lower_weight + b.z * upper_weight,
			a.w * lower_weight + b.w * upper_weight);
	}
	else
	{
		return a;
	}
}

#pragma mark - ---

inline vec2::vec2(vec3 const& v) : x(v.x), y(v.y) {}

inline vec3::vec3(vec4 const& v) : x(v.x), y(v.y), z(v.z) {}

inline mat3::mat3(mat4 const& m) : x(m.x), y(m.y), z(m.z) {}

inline mat3::mat3(quat const& q)
{
	float wx, wy, wz, xx, yy, yz, xy, xz, zz, x2, y2, z2;
	float s = 2 / q.norm();

	x2 = q.x * s;   y2 = q.y * s;   z2 = q.z * s;
	xx = q.x * x2;  xy = q.x * y2;  xz = q.x * z2;
	yy = q.y * y2;  yz = q.y * z2;  zz = q.z * z2;
	wx = q.w * x2;  wy = q.w * y2;  wz = q.w * z2;

	x.x = 1 - (yy + zz);
	y.x = xy - wz;
	z.x = xz + wy;

	x.y = xy + wz;
	y.y = 1 - (xx + zz);
	z.y = yz - wx;

	x.z = xz - wy;
	y.z = yz + wx;
	z.z = 1 - (xx + yy);
}
