from gaussian_gradient_magnitude import run as run_sanders_ggm

location = '/Users/emoravec/Documents/Research/merging_clusters/analysis/MOO_1142/'

mooNet = location+'xray/images/XMM/XMM_comb-net-image.fits'
mooNetCenter = location+'xray/images/XMM/XMM_comb-net-center.fits' 
mooObj = location+'xray/images/XMM/XMM_comb-obj-im-400-7200.fits'
for kernel in range(1, 6):
    print(f"Processing scale {kernel}")
    run_sanders_ggm(infile=mooNetCenter,
                    outfile=location+'ggm/XMM/sanders/moo1142_xxm.combNetCenter.sanders.ggm.{0}.fits'.format(kernel),
                    scale=kernel)
    # run_sanders_ggm(infile=mooObj,
    #                 outfile=location+'ggm/XMM/sanders/moo1142_xxm.combObj.sanders.ggm.{0}.fits'.format(kernel),
    #                 scale=kernel)