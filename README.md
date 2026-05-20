# Omicslog

## Importing packages


```python
from omicslog import log_start
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix

```

<div style="background-color:#e7f3fe; border-left:6px solid #2196F3; padding:12px 16px; border-radius:4px; margin:8px 0;">
<strong>📝 Note</strong><br>
The AnnData object were generated using code from the original <a href="https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html"> AnnData </a> documentation.
</div>


```python
counts = csr_matrix(np.random.poisson(1, size=(100, 2000)), dtype=np.float32)
adata = ad.AnnData(counts)

adata.obs_names = [f"Cell_{i:d}" for i in range(adata.n_obs)]
adata.var_names = [f"Gene_{i:d}" for i in range(adata.n_vars)]

logdata = log_start(adata)
print(logdata)

ct = np.random.choice(["B", "T", "Monocyte"], size=(logdata.n_obs,))
logdata.obs["cell_type"] = pd.Categorical(ct)  # Categoricals are preferred for efficiency
print(logdata)
logdata.uns["_omicslog"]
```

    AnnData object with n_obs × n_vars = 100 × 2000
        uns: '_omicslog'
    AnnData object with n_obs × n_vars = 100 × 2000
        obs: 'cell_type'
        uns: '_omicslog'
    
    Operation log:
    [2026-05-20 13:49:41] obs: 'cell_type' added





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Time</th>
      <th>Operation</th>
      <th>Message</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2026-05-20 13:49:41</td>
      <td>obs</td>
      <td>'cell_type' added</td>
    </tr>
  </tbody>
</table>
</div>



## Fltrating by Cells (.obs)


```python
logdata = logdata[logdata.obs.cell_type == "B"]
print(logdata)
logdata.uns["_omicslog"]
```

    AnnData object with n_obs × n_vars = 27 × 2000
        obs: 'cell_type'
        uns: '_omicslog'
    
    Operation log:
    [2026-05-20 13:49:41] obs: 'cell_type' added
    [2026-05-20 13:49:43] subset: removed 73 samples (73%), 27 samples remaining





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Time</th>
      <th>Operation</th>
      <th>Message</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2026-05-20 13:49:41</td>
      <td>obs</td>
      <td>'cell_type' added</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2026-05-20 13:49:43</td>
      <td>subset</td>
      <td>removed 73 samples (73%), 27 samples remaining</td>
    </tr>
  </tbody>
</table>
</div>



## Filtering by Genes (.var)


```python
logdata = logdata[:,logdata.var_names.str.endswith("1")]
print(logdata)
logdata.uns["_omicslog"]
```

    AnnData object with n_obs × n_vars = 27 × 200
        obs: 'cell_type'
        uns: '_omicslog'
    
    Operation log:
    [2026-05-20 13:49:41] obs: 'cell_type' added
    [2026-05-20 13:49:43] subset: removed 73 samples (73%), 27 samples remaining
    [2026-05-20 13:49:46] subset: removed 1800 genes (90%), 200 genes remaining





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Time</th>
      <th>Operation</th>
      <th>Message</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2026-05-20 13:49:41</td>
      <td>obs</td>
      <td>'cell_type' added</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2026-05-20 13:49:43</td>
      <td>subset</td>
      <td>removed 73 samples (73%), 27 samples remaining</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2026-05-20 13:49:46</td>
      <td>subset</td>
      <td>removed 1800 genes (90%), 200 genes remaining</td>
    </tr>
  </tbody>
</table>
</div>



## Adding observatons and variables


```python
logdata.obsm["X_umap"] = np.random.normal(0, 1, size=(logdata.n_obs, 2))
logdata.varm["gene_stuff"] = np.random.normal(0, 1, size=(logdata.n_vars, 5))
print(logdata)
logdata.uns["_omicslog"]
```

    AnnData object with n_obs × n_vars = 27 × 200
        obs: 'cell_type'
        uns: '_omicslog'
        obsm: 'X_umap'
        varm: 'gene_stuff'
    
    Operation log:
    [2026-05-20 13:49:41] obs: 'cell_type' added
    [2026-05-20 13:49:43] subset: removed 73 samples (73%), 27 samples remaining
    [2026-05-20 13:49:46] subset: removed 1800 genes (90%), 200 genes remaining
    [2026-05-20 13:49:48] obsm: 'X_umap' added
    [2026-05-20 13:49:48] varm: 'gene_stuff' added





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Time</th>
      <th>Operation</th>
      <th>Message</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2026-05-20 13:49:41</td>
      <td>obs</td>
      <td>'cell_type' added</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2026-05-20 13:49:43</td>
      <td>subset</td>
      <td>removed 73 samples (73%), 27 samples remaining</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2026-05-20 13:49:46</td>
      <td>subset</td>
      <td>removed 1800 genes (90%), 200 genes remaining</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2026-05-20 13:49:48</td>
      <td>obsm</td>
      <td>'X_umap' added</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2026-05-20 13:49:48</td>
      <td>varm</td>
      <td>'gene_stuff' added</td>
    </tr>
  </tbody>
</table>
</div>



## Adding layers


```python
logdata.layers["log_transformed"] = np.log1p(logdata.X)
print(logdata)
logdata.uns["_omicslog"]
```

    AnnData object with n_obs × n_vars = 27 × 200
        obs: 'cell_type'
        uns: '_omicslog'
        obsm: 'X_umap'
        varm: 'gene_stuff'
        layers: 'log_transformed'
    
    Operation log:
    [2026-05-20 13:49:41] obs: 'cell_type' added
    [2026-05-20 13:49:43] subset: removed 73 samples (73%), 27 samples remaining
    [2026-05-20 13:49:46] subset: removed 1800 genes (90%), 200 genes remaining
    [2026-05-20 13:49:48] obsm: 'X_umap' added
    [2026-05-20 13:49:48] varm: 'gene_stuff' added
    [2026-05-20 13:49:50] layers: 'log_transformed' added





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Time</th>
      <th>Operation</th>
      <th>Message</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>2026-05-20 13:49:41</td>
      <td>obs</td>
      <td>'cell_type' added</td>
    </tr>
    <tr>
      <th>1</th>
      <td>2026-05-20 13:49:43</td>
      <td>subset</td>
      <td>removed 73 samples (73%), 27 samples remaining</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2026-05-20 13:49:46</td>
      <td>subset</td>
      <td>removed 1800 genes (90%), 200 genes remaining</td>
    </tr>
    <tr>
      <th>3</th>
      <td>2026-05-20 13:49:48</td>
      <td>obsm</td>
      <td>'X_umap' added</td>
    </tr>
    <tr>
      <th>4</th>
      <td>2026-05-20 13:49:48</td>
      <td>varm</td>
      <td>'gene_stuff' added</td>
    </tr>
    <tr>
      <th>5</th>
      <td>2026-05-20 13:49:50</td>
      <td>layers</td>
      <td>'log_transformed' added</td>
    </tr>
  </tbody>
</table>
</div>


