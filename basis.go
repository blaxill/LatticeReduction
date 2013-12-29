package LatticeReduction

import (
	"errors"
	// "fmt"
	// "math/big"
)

var (
	ErrInconsistantDimension  = errors.New("Basis vectors have inconsistant dimension.")
	ErrNotLinearlyIndependant = errors.New("Basis vectors are not linearly independant.")
)


type Basis interface{
	Rank() int
	Dimension() int
	Copy() Basis

	// Dot(column1, column2 int, out *big.Int)
	FDot(column1, column2 int)float64
	// ColumnReduce(column1, column2 int, mu *big.Int) // 
	FColumnReduce(column1, column2 int, mu float64) // column1 -= mu*column2

	ColumnSwap(column1,column2 int)

	// Get(column,row int, out *big.Int)
	FGet(column,row int)float64
}

func dot(lhs , rhs []float64) (r float64) {
	for i, x := range lhs {
		r += x * rhs[i]
	}
	return
}

func abs(v float64) float64 {
	if v < 0 {
		v = -v
	}

	return v
}