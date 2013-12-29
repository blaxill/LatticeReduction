package LatticeReduction

// Algorithm from -
// Lattice Basis Reduction:
// Improved Practical Algorithms and Solving Subset Sum Problems

import (
	"math"
)

// L3FP is L^3-Reduction in floating point.
// Basis is kept in exact representation, but orthogonal basis
// and friends are in floating point.
// Delta should be 0.5 < delta < 1
func L3FP(inBasis Basis, delta float64) Basis {
	var (
		basis = inBasis.Copy()
		k     = 1
		Fc    = false
		bd    = make([][]float64, basis.Rank())
		mu    = make([][]float64, basis.Rank())
		c     = make([]float64, basis.Rank())
	)

	for i := 0; i<basis.Rank();i++ {
		bd[i] = make([]float64, basis.Dimension())
		mu[i] = make([]float64, basis.Rank())
		for j := 0; j < basis.Dimension(); j++ {
			bd[i][j] = basis.FGet(i,j)
		}
	}

	for k < basis.Rank() {
		// 2.
		c[k] = dot(bd[k], bd[k])
		if k == 1 {
			c[0] = dot(bd[0], bd[0])
		}
		for j := 0; j < k; j++ {
			if abs(dot(bd[k], bd[j])) < _2_p_nh_tor*math.Sqrt(dot(bd[k], bd[k])*dot(bd[j], bd[j])) {
				mu[k][j] = basis.FDot(k,j)
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
				basis.FColumnReduce(k, j,_mu)
				for i := range bd[k] {
					bd[k][i] = basis.FGet(k,i)
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

		// 4.
		if delta*c[k-1] > c[k]+mu[k][k-1]*mu[k][k-1]*c[k-1] {

			basis.ColumnSwap(k,k-1)
			bd[k], bd[k-1] = bd[k-1], bd[k]

			k -= 1
			if k < 1 {
				k = 1
			}
		} else {
			k += 1
		}
	}
	return basis
}