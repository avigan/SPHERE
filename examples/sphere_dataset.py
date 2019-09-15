import vltpf.SPHERE as SPHERE

#%% init data set
ds = SPHERE.Dataset('/Users/avigan/data/VLTPF-test-target/', log_level='info')

print('IRDIS reductions:')
for red in ds.IRDIS_reductions:
    print(red)
    red.config['clean'] = True
print()

print('IFS reductions:')
for red in ds.IFS_reductions:
    print(red)
    red.config['clean'] = True
print()

#%% full reduction with default parameters
ds.full_reduction()
