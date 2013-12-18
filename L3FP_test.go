package LatticeReduction

import (
	"fmt"
	"math/rand"
	"testing"
)

var SmallBasisTest = Basis{
	Int64Vector{1, 1, 1},
	Int64Vector{-1, 0, 2},
	Int64Vector{3, 5, 6},
}

func BenchmarkSmallBasisL3FP100x100(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		basis := make(Basis, 100)
		for j := 0; j < 100; j++ {
			v := make(Int64Vector, 100)
			basis[j] = v
			for k := 0; k < 100; k++ {
				v[k] = rand.Int63n(0x7FFFFFFFFFFFFFF)
			}
		}
		b.StartTimer()
		_ = basis.L3FP(0.99)
	}
}

func ExampleReduceL3FP() {
	fmt.Println(SmallBasisTest)
	fmt.Println(SmallBasisTest.L3FP(0.9))
	// Output:
	// [1 1 1]
	// [-1 0 2]
	// [3 5 6]
	//
	// [0 1 0]
	// [1 0 1]
	// [-1 0 2]
}

