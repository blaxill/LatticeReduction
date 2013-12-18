package LatticeReduction

import (
	"errors"
	"fmt"
)

var (
	ErrInconsistantDimension  = errors.New("Basis vectors have inconsistant dimension.")
	ErrNotLinearlyIndependant = errors.New("Basis vectors are not linearly independant.")
)

type Float64Vector []float64

type Vector interface {
	// FDot(Float64Vector)float64
	Dot(Vector)float64
	FReduceInplace(float64,Vector)
	Len() int
	Copy() Vector
	FAt(int) float64
}

type Basis []Vector

func (lhs Float64Vector) FDot(rhs Float64Vector)(r float64) {
	for i, x := range lhs {
		r += x * rhs[i]
	}
	return
}

func dot(lhs Float64Vector, rhs Float64Vector)(r float64){
			for i, x := range lhs {
				r += x * rhs[i]
			}
			return
}

func abs(v float64)float64{
			if v < 0 {
				v = -v
			}

			return v
}

func (b Basis) Copy() Basis {
	o := make([]Vector, len(b))
	for i, v := range b {
		o[i] = v.Copy()
	}
	return o
}

func (b Basis) String() (s string) {
	for _, r := range b {
		s += fmt.Sprintln(r)
	}
	return s
}
