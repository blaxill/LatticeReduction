package LatticeReduction

import (
	"fmt"
	"math/big"
)

// Lattice basis for int64 sized values
type Int64Basis [][]int64

func (b Int64Basis) Copy() Basis {
	o := make(Int64Basis, len(b))
	for i, x := range b {
		o[i] = make([]int64, len(x))
		for j, y := range x {
			o[i][j] = y
		}
	}
	return o
}

func (b Int64Basis) ColumnSwap(column1, column2 int) {
	b[column1], b[column2] = b[column2], b[column1]
}

func (b Int64Basis) Rank() int {
	return len(b)
}

func (b Int64Basis) Dimension() int {
	return len(b[0])
}

func (b Int64Basis) FDot(column1, column2 int) (r float64) {
	for i, x := range b[column1] {
		r += float64(x) * float64(b[column2][i])
	}
	return
}

func (b Int64Basis) FPairSize(column1, column2 int) (sub,add float64) {
	for i, x := range b[column1] {
		ts:=float64(x - b[column2][i])
		ta:=float64(x + b[column2][i])
		sub += ts*ts
		add += ta*ta 
	}
	return
}

func (b Int64Basis) ColumnReduceInt64(column1, column2 int, mu int64) {
	for i := range b[column1] {
		b[column1][i] -= mu * b[column2][i]
	}
}

func (b Int64Basis) FGet(column, row int) float64 {
	return float64(b[column][row])
}

func (b Int64Basis) String() (s string) {
	for _, r := range b {
		s += fmt.Sprintln(r)
	}
	return s
}

func (b Int64Basis) PremoteToBig() Basis {
	o := make(BigBasis, len(b))
	for i, x := range b {
		o[i] = make([]*big.Int, len(x))
		for j, y := range x {
			o[i][j] = big.NewInt(y)
		}
	}
	return o
}

func (b Int64Basis) Mod() (r float64) {

	for _, x := range b {
		for _, y := range x {
			r += float64(y) * float64(y)
		}
	}
	r /= float64(b.Rank()) * float64(b.Dimension())
	return
}
