import plur_cert_positive as certificate
certificate.RHO_MIN = 0.4
certificate.RHO_MAX = 0.42
raise SystemExit(certificate.main([
    "--rho-min", "0.4", "--rho-max", "0.42",
    "--rho-panels", "4", "--x-panels", "3200",
    "--min-rho-width", "0.000025"
]))
