import vltpf.SPHERE as SPHERE

#%% init data set
ds = SPHERE.Dataset('/Users/avigan/data/VLTPF-test-target/test/', log_level='info')

print('IRDIS reductions:')
for red in ds.IRDIS_reductions:
    print(red)
print()

print('IFS reductions:')
for red in ds.IFS_reductions:
    print(red)
print()

#%% full reduction with default parameters
ds.full_reduction()
