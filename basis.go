package LatticeReduction

import "errors"

// SmallBasis is for when all elements of a basis can be represented by int64 types.
// For stability reasons, elements should be less than 2^59
type SmallBasis [][]int64

var (
	ErrInconsistantDimension  = errors.New("Basis vectors have inconsistant dimension.")
	ErrNotLinearlyIndependant = errors.New("Basis vectors are not linearly independant.")
)

// Validate
// TODO: check for linear dependance
func (b SmallBasis) Validate() error {
	var m = len(b)
	if m == 0 {
		return nil
	}

	var n = len(b[0])

	if n < m {
		return ErrNotLinearlyIndependant
	}

	for i := range b {
		if len(b[i]) != n {
			return ErrInconsistantDimension
		}
	}

	return nil
}

func (b SmallBasis) Copy() SmallBasis {
	o := make([][]int64, len(b))
	for i, v := range b {
		o[i] = make([]int64, len(b[i]))
		copy(o[i], v)
	}
	return o
}
