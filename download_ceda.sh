#
# Download ScenarioMIP (CMIP6) data for selected models and scenarios
# from the CEDA data catalogue at https://data.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP
#

for model in "CAS/FGOALS-g3" "MIROC/MIROC6" "NCAR/CESM2-WACCM" "AWI/AWI-CM-1-1-MR"; do
    for ssp in "ssp126" "ssp245" "ssp370" "ssp585"; do
        for variable in "tas" "pr"; do
            echo "SSP: $ssp  VAR: $variable"
            echo "wget -e robots=off --mirror --no-parent --no-clobber -r https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/$model/$ssp/r1i1p1f1/day/$variable/gn/latest/"
            wget -e robots=off --no-parent --no-clobber -r "https://dap.ceda.ac.uk/badc/cmip6/data/CMIP6/ScenarioMIP/$model/$ssp/r1i1p1f1/day/$variable/gn/latest/"
        done
    done
done
