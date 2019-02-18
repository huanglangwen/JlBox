
function lognormal(num_bins::Integer,total_conc::Real,meansize::Real,
                   size_std::Real,lowersize::Real,uppersize::Real)
    mu=log(meansize)
    sigma=log(size_std)
    x_output=exp.(linspace(log(lowersize),log(uppersize),num_bins))
    pdf=x->exp(-(log(x)-mu)^2/(2*sigma^2))/(x*sigma*sqrt(2*pi))
    pdf_output=[pdf(x) for x in x_output]
    return (pdf_output./sum(pdf_output))*total_conc,x_output
end
