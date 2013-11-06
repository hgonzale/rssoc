function [ c_time, cost ] = cost_and_computation_time( result )

c_time = 0;

for i = 1:length( result.optfctn )
    if ( ~isempty( result.optfctn( i ).time ) )
        c_time = c_time + result.optfctn( i ).time;
    end
end
for i = 1:length( result.res_insmode )
    if ( ~isempty( result.res_insmode( i ).time ) )
        c_time = c_time + result.res_insmode( i ).time;
    end
end
for i = 1:length( result.res_pwm )
    if ( ~isempty( result.res_pwm( i ).time ) )
        c_time = c_time + result.res_pwm( i ).time;
    end
end

index = length( result.res_pwm );

cost = obj_fctn( result.user( index ),  result.res_pwm( index ).tau, result.res_pwm( index ).u, result.res_pwm( index ).d );