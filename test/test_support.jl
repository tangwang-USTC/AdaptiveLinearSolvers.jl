function standard_matrix()
    [4.0 1.0; 1.0 3.0]
end

standard_rhs() = [1.0, 2.0]

function standard_spd_contract()
    MathematicalContract(
        square=PropertyEvidence(Certified; source=:caller),
        hermitian=PropertyEvidence(Certified; source=:caller),
        positive_definite=PropertyEvidence(Certified; source=:caller),
    )
end
