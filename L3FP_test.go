package LatticeReduction

import (
	"fmt"
	"math/rand"
	"testing"
)

var (
	SmallBasisTest = Int64Basis{
		[]int64{1, 1, 1},
		[]int64{-1, 0, 2},
		[]int64{3, 5, 6},
	}

	LargeBasisTest = SmallBasisTest.PremoteToBig()
)

const size = 35

func TestXx(t *testing.T) {
	cases:=250
	if testing.Short(){
		cases=10
	}
	count:=0
	for i := 0; i < cases; i++ {
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0xFFFFFF)
			}
		}
		notDeep,deep:=L3FP(basis, 0.75).(Int64Basis).Mod(), L3FPDeep(basis, 0.75,4).(Int64Basis).Mod()
		if notDeep<=deep{
			count++
		}
	}
	t.Logf("L3FP %d: L3FPDeep %d  (events smaller)",count,cases-count)
}


func TestXxx(t *testing.T) {
	cases:=250
	if testing.Short(){
		cases=10
	}
	count:=0
	for i := 0; i < cases; i++ {
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0xFFFFFF)
			}
		}
		notDeep,deep:=L3FP(basis, 0.75).(Int64Basis).Mod(), PairwiseReduce(L3FP(basis, 0.75)).(Int64Basis).Mod()
		if notDeep<=deep{
			count++
		}
	}
	t.Logf("L3FP %d: Pairwise L3FP %d  (events smaller)",count,cases-count)
}


func TestXxxx(t *testing.T) {
	cases:=250
	if testing.Short(){
		cases=10
	}
	count:=0
	for i := 0; i < cases; i++ {
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0xFFFFFF)
			}
		}
		notDeep,deep:=L3FPDeep(basis, 0.75,4).(Int64Basis).Mod(), PairwiseReduce(L3FPDeep(basis, 0.75,4)).(Int64Basis).Mod()
		if notDeep<=deep{
			count++
		}
	}
	t.Logf("L3FPDeep %d: Pairwise L3FPDeep %d  (events smaller)",count,cases-count)
}

func BenchmarkSmallBasisL3FP20x20(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0x7FFFFFFFFFFFFFF)
			}
		}
		b.StartTimer()
		_ = L3FP(basis, 0.75)
	}
}

func BenchmarkSmallBasisL3FPDeep20x20(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0x7FFFFFFFFFFFFFF)
			}
		}
		b.StartTimer()
		_ = L3FPDeep(basis, 0.75, 4)
	}
}

func BenchmarkLargeBasisL3FP20x20(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0x7FFFFFFFFFFFFFF)
			}
		}
		large := basis.PremoteToBig()
		b.StartTimer()
		_ = L3FP(large, 0.75)
	}
}

func BenchmarkLargeBasisL3FPDeep20x20(b *testing.B) {
	for i := 0; i < b.N; i++ {
		b.StopTimer()
		basis := make(Int64Basis, size)
		for j := 0; j < size; j++ {
			v := make([]int64, size)
			basis[j] = v
			for k := 0; k < size; k++ {
				v[k] = rand.Int63n(0x7FFFFFFFFFFFFFF)
			}
		}
		large := basis.PremoteToBig()
		b.StartTimer()
		_ = L3FPDeep(large, 0.75,4)
	}
}

func ExampleReduceSmallL3FP() {
	fmt.Println(SmallBasisTest)
	fmt.Println(L3FP(SmallBasisTest, 0.75))
	// Output:
	// [1 1 1]
	// [-1 0 2]
	// [3 5 6]
	//
	// [0 1 0]
	// [1 0 1]
	// [-1 0 2]
}

func ExampleReduceBigL3FP() {
	fmt.Println(LargeBasisTest)
	fmt.Println(L3FP(LargeBasisTest, 0.75))
	// Output:
	// [1 1 1]
	// [-1 0 2]
	// [3 5 6]
	//
	// [0 1 0]
	// [1 0 1]
	// [-1 0 2]
}

func ExampleReduceSmallL3FPDeep() {
	fmt.Println(SmallBasisTest)
	fmt.Println(L3FPDeep(SmallBasisTest, 0.75, 4))
	// Output:
	// [1 1 1]
	// [-1 0 2]
	// [3 5 6]
	//
	// [0 1 0]
	// [1 0 1]
	// [-1 0 2]
}

func ExampleReduceBigL3FPDeep() {
	fmt.Println(LargeBasisTest)
	fmt.Println(L3FPDeep(LargeBasisTest, 0.75, 4))
	// Output:
	// [1 1 1]
	// [-1 0 2]
	// [3 5 6]
	//
	// [0 1 0]
	// [1 0 1]
	// [-1 0 2]
}
