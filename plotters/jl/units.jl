const TIME_STEP_AU::Float64 = 1. * 1e10 
const STEPS_PER_SAVE::Int32 = 1000

function time_au_to_us(time_au::Float64)
    return 2.418884 * 1e-11 * time_au
end

function time_au_to_ms(time_au::Float64)
    return 2.418884 * 1e-14 * time_au
end
