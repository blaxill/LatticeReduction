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

type Basis interface {
	Rank() int
	Dimension() int
	Copy() Basis

	// TODO: Dot(column1, column2 int, out *big.Int)

	// column1*column2 as float64
	FDot(column1, column2 int) float64
	// Dot(column1-column2), Dot(column1+column2)
	FPairSize(column1, column2 int) (float64,float64)

	// TODO: ColumnReduce(column1, column2 int, mu *big.Int)

	// column1 -= mu*column2
	ColumnReduceInt64(column1, column2 int, mu int64)



	ColumnSwap(column1, column2 int)

	// TODO: Get(column,row int, out *big.Int)

	FGet(column, row int) float64
}

func dot(lhs, rhs []float64) (r float64) {
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
