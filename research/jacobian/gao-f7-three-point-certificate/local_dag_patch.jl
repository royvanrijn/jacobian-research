# Adapted from Groebner.jl 0.10.6 src/f4/matrix.jl, GNU GPL v2.
# Original source hash and package lock are retained in dag-provenance.json.
# Only the global expansion is replaced: emit local identities to older basis
# elements, keep empty global change entries, stop immediately on a constant.
const F7_DAG_IO=Ref{IO}(stdout)
const F7_DAG_INITIAL=Ref(false)
function f7_emit_node(idx,basis,ht,relations)
    io=F7_DAG_IO[]
    println(io,"NODE ",idx)
    for j in eachindex(basis.monoms[idx])
        e=zeros(Int,ht.nvars);monom_to_vector!(e,ht.monoms[basis.monoms[idx][j]])
        println(io,"P ",basis.coeffs[idx][j],"|",join(e,","))
    end
    if !isnothing(relations)
        values=Tuple{Int,Tuple,Any}[]
        for ((parent,monom),cf) in relations
            iszero(cf) && continue
            e=zeros(Int,ht.nvars);monom_to_vector!(e,ht.monoms[monom])
            push!(values,(parent,Tuple(e),cf))
        end
        sort!(values,by=x->(x[1],x[2]))
        for (parent,e,cf) in values
            @assert parent<idx
            println(io,"R ",parent,"|",cf,"|",join(e,","))
        end
    end
    println(io,"END")
end
function matrix_convert_rows_to_basis_elements!(
    matrix::MacaulayMatrix,
    basis::Basis{C},
    ht::MonomialHashtable,
    symbol_ht::MonomialHashtable,
    params::AlgorithmParameters;
    batched_ht_insert::Bool=false
) where {C}
    if !F7_DAG_INITIAL[]
        for idx in 1:basis.n_processed
            f7_emit_node(idx,basis,ht,nothing)
        end
        F7_DAG_INITIAL[]=true
    end
    # We mutate the basis directly by adding new elements

    basis_resize_if_needed!(basis, matrix.npivots)
    rows = matrix.lower_rows
    crs = basis.n_processed

    _, _, nl, nr = matrix_block_sizes(matrix)
    support_size = 0
    for i in 1:length(rows)
        if isassigned(rows, i)
            support_size += length(rows[i])
        end
    end

    if batched_ht_insert && matrix.npivots > 1 && support_size > 4 * nr
        column_to_basis_monom = zeros(MonomId, nr)
        @inbounds for i in 1:(matrix.npivots)
            for j in 1:length(rows[i])
                idx = rows[i][j]
                column_to_basis_monom[idx - nl] = 1
            end
        end
        matrix_insert_in_basis_hashtable_pivots_masked!(
            column_to_basis_monom,
            ht,
            symbol_ht,
            matrix.column_to_monom,
            nl
        )
        @inbounds for i in 1:(matrix.npivots)
            colidx = rows[i][1]
            for j in 1:length(rows[i])
                columnid = rows[i][j]
                rows[i][j] = column_to_basis_monom[columnid - nl]
            end
            basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[colidx]]
            basis.monoms[crs + i] = matrix.lower_rows[i]
            @invariant length(basis.coeffs[crs + i]) == length(basis.monoms[crs + i])
        end
    else
        @inbounds for i in 1:(matrix.npivots)
            colidx = rows[i][1]
            matrix_insert_in_basis_hashtable_pivots!(rows[i], ht, symbol_ht, matrix.column_to_monom)
            basis.coeffs[crs + i] = matrix.some_coeffs[matrix.lower_to_coeffs[colidx]]
            basis.monoms[crs + i] = matrix.lower_rows[i]
            @invariant length(basis.coeffs[crs + i]) == length(basis.monoms[crs + i])
        end
    end

    if params.changematrix
        resize!(basis.changematrix, basis.n_filled + matrix.npivots)
        for i in 1:matrix.npivots
            idx=crs+i
            basis.changematrix[idx]=Dict{Int,Dict{MonomId,C}}()
            f7_emit_node(idx,basis,ht,matrix.changematrix[i])
            if length(basis.monoms[idx])==1 && monom_totaldeg(ht.monoms[basis.monoms[idx][1]])==0
                flush(F7_DAG_IO[])
                error("F7_UNIT_EXPORTED")
            end
        end
    end

    basis.n_filled += matrix.npivots
end

