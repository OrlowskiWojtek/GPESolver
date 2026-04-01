const TIME_STEP_AU::Float64 = 1. * 1e10 
const STEPS_PER_SAVE::Int32 = 1000

function time_au_to_us(time_au::Float64)
    return 2.418884 * 1e-11 * time_au
end

function time_au_to_ms(time_au::Float64)
    return 2.418884 * 1e-14 * time_au
end

function length_au_to_nm(length_au::Float64)
    return length_au * 0.052917721092
end

function length_au_to_μm(length_au::Float64)
    return length_au * 0.000052917721092
end

function length_nm_to_au(length_nm::Float64)
    return length_nm / 0.052917721092
end

function length_nm_to_μm(length_nm::Float64)
    return length_nm * 0.001
end

function length_nm3_to_μm3(length_nm::Float64)
    return length_nm * 0.001^3
end

function length_au2_to_μm2(length_au::Float64)
    return length_au * 0.000052917721092^2
end

function length_au3_to_μm3(length_au::Float64)
    return length_au * 0.000052917721092^3
end
