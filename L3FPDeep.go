package LatticeReduction

// Algorithm from -
// Lattice Basis Reduction:
// Improved Practical Algorithms and Solving Subset Sum Problems

import (
	"math"
)

// L3FPDeep is L3FP with 'deep insertion'
// This is super-polynomial time in worst case, but will output a more reduced basis.
// Delta should be 0.5 < delta < 1
func (inBasis SmallBasis) L3FPDeep(delta float64) SmallBasis {
	var (
		basis = inBasis.Copy()
		k     = 1
		Fc    = false
		bd    = make([][]float64, len(basis))
		mu    = make([][]float64, len(basis))
		c     = make([]float64, len(basis))

		dot = func(lhs []float64, rhs []float64) (r float64) {
			for i, x := range lhs {
				r += x * rhs[i]
			}
			return
		}
		idot = func(lhs []int64, rhs []int64) (r float64) {
			for i, x := range lhs {
				r += float64(x * rhs[i])
			}
			return
		}

		abs = func(v float64) float64 {
			if v < 0 {
				v = -v
			}

			return v
		}

		vecReduce = func(lhs []int64, mu float64, rhs []int64) {
			for i := range lhs {
				lhs[i] -= int64(mu * float64(rhs[i]))
			}
		}
	)

	for i, v := range basis {
		bd[i] = make([]float64, len(v))
		mu[i] = make([]float64, len(basis))
		for j := range v {
			bd[i][j] = float64(v[j])
		}
	}

L:
	for k < len(basis) {
		// 2.
		c[k] = dot(bd[k], bd[k])
		if k == 1 {
			c[0] = dot(bd[0], bd[0])
		}
		for j := 0; j < k; j++ {
			if abs(dot(bd[k], bd[j])) < _2_p_h_tor*math.Sqrt(dot(bd[k], bd[k])*dot(bd[j], bd[j])) {
				mu[k][j] = idot(basis[k], basis[j])
			} else {
				mu[k][j] = dot(bd[k], bd[j])
			}

			for i := 0; i < j; i++ {
				mu[k][j] -= mu[j][i] * mu[k][i] * c[i]
			}

			mu[k][j] /= c[j]

			c[k] -= mu[k][j] * mu[k][j] * c[j]
		}

		// 3.
		for j := k - 1; j >= 0; j-- {
			if abs(mu[k][j]) > 0.5 {
				var _mu = math.Floor(mu[k][j] + 0.5)

				if abs(_mu) > _2_p_h_tor {
					Fc = true
				}

				for i := 0; i < j; i++ {
					mu[k][i] -= _mu * mu[j][i]
				}

				mu[k][j] -= _mu
				vecReduce(basis[k], _mu, basis[j])
				for i := range bd[k] {
					bd[k][i] = float64(basis[k][i])
				}
			}
		}
		if Fc {
			Fc = false
			k -= 1
			if k < 1 {
				k = 1
			}
			continue
		}

		// New step 4.
		_c := dot(bd[k], bd[k])
		for i := 0; i < k; i++ {
			if delta*c[i] <= _c {
				_c -= mu[k][i] * mu[k][i] * c[i]
			} else {
				s := basis[len(basis)-1]
				basis = append(basis, nil)
				copy(basis[i+1:len(basis)-1], basis[i:len(basis)-2])
				basis[i] = s
				basis = basis[:len(basis)-1]

				sd := bd[len(bd)-1]
				bd = append(bd, nil)
				copy(bd[i+1:len(bd)-1], bd[i:len(bd)-2])
				bd[i] = sd
				bd = bd[:len(bd)-1]

				k -= 1
				if k < 1 {
					k = 1
				}
				continue L
			}
		}
		k += 1
	}
	return basis
}
