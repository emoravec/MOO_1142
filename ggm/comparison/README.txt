Write up of tests

in Ds9 XMM_comb-net-center.fits (from Craig) = XMM_comb-obj-im-400-7200.fits (from Tony) -> cut out of cluster
XMM_comb-net-image.fits (from Tony) -> the whole field
See XMM_orig_compare.png

Scipy results
XMM_comb-net-center.fits & XMM_comb-net-image.fits => the point sources are large zero points
XMM_comb-obj-im-400-7200.fits -> don't have zero problem but PN chip lines are seen

Sanders results
Same -> see compare_sanders_combObj_center.png -> left: XMM_comb-obj-im-400-7200.fits and right: XMM_comb-obj-im-400-7200.fits

Makes sense - it is using the same underlying function scipy.ndimage.gaussian_gradient_magnitude(img, scale)

Tony gave me XMM_comb-net-image.fits which is the whole field and XMM_comb-obj-im-400-7200.fits which is a cutout version of that. When I run XMM_comb-net-image.fits through either scipy.ndimage.gaussian_gradient_magnitude(img, scale) OR using jeremey's code, the areas that are zeros due to point source subtraction become larger and larger sections of zeros as the kernel size increases (see moo_1142_xxm.combNetCenter.ggm.all.pdf and moo_1142_xxm.combNetObj.ggm.all.pdf). But XMM_comb-obj-im-400-7200.fits does not have this issue (see moo_1142_xxm.combNetCenter.ggm.all.pdf). However, the chip extents are quite obvious. 

I thought to myself, Craig's image didn't look like that (see his post on December 27th and the attached xmm_ggm_2_sarazin.png). So today I talked to Craig and he gave me the image he was using XMM_comb-net-center.fits. So I ran the image he was using through both scipy.ndimage.gaussian_gradient_magnitude(img, scale) and jeremey's code get the same issue with the increasing size of zeros as the kernel size increases (see compare_sanders_combObj_center.png -> left: XMM_comb-net-center.fits from Craig and right: XMM_comb-obj-im-400-7200.fits from Tony). 

I have no idea what is going on here and not sure how to procede. The only thing I can think of is different versions of scipy?? But that would be quite bad if the difference was so big about how the function behaved between scipy version.

Anyway, so right now I cannot reproduce the image Craig made and am going to work with the ggm he sent me (XMM_ggm_2_sarazin.fits)and the one I've made that has some holes to compare the M2, Chandra, and XMM ggms and try to pin down the location of the shock.