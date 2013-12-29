package LatticeReduction

func PairwiseReduce(inBasis Basis) Basis{
	var (
		basis = inBasis.Copy()
		)

	L:
	for i := 0; i < basis.Rank(); i++{
		v:=basis.FDot(i,i)
		for k := 0; k < basis.Rank(); k++{
			if k == i {
				continue
			}

			s,a:=basis.FPairSize(i,k)
			m:=1
			if a<s{
				m=-1
				s=a
			}

			if s<v{
				basis.ColumnReduceInt64(i,k,int64(m))

				goto L
			}
		}
	}

	return basis
}