
## Use Cases
Based on your study you might be interested in one of the following use cases.

### UseCase 1: qPCR experiment with multiple reference genes ###
You have multiple reference genes (usually 2 or 3) and their corresponding raw CT values. Here is how you can optimally aggregate them into one new internal control using weighted geometric mean:

**step 1)** The CT values can be provided in one of these two formats:  

* a tab separated .txt file:

```
      Sample1  Sample2 Sample3 Sample4 Sample5 Sample6 Sample7 Sample8 Sample9 Sample10 Sample11 Sample12
RNU44 25.8800  24.338  24.666  24.224 23.7500  24.642  23.720 23.2620  24.318   23.624   24.930   24.104
RNU48 21.3680  21.756  23.174  22.252 21.7820  22.994  21.564 20.9600  21.886   22.890   23.618   22.060
RNU6B 32.2575  29.480  30.455  29.765 29.7175  30.015  29.735 29.5475  30.640   30.545   30.755   30.335
```
* An excel file:

 \s | Sample1 | Sample2 | Sample3 | Sample4 | Sample5 | Sample6 | Sample7 | Sample8 | Sample9 | Sample10 | Sample11 | Sample12
------| -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | --------
RNU44 | 25.88 | 24.338 | 24.666 | 24.224 | 23.75 | 24.642 | 23.72 | 23.262 | 24.318 | 23.624 | 24.93 | 24.104
RNU48 | 21.368 | 21.756 | 23.174 | 22.252 | 21.782 | 22.994 | 21.564 | 20.96 | 21.886 | 22.89 | 23.618 | 22.06
RNU6B | 32.2575 | 29.48 | 30.455 | 29.765 | 29.7175 | 30.015 | 29.735 | 29.5475 | 30.64 | 30.545 | 30.755 | 30.335

**step 2)**
If your data only contains reference genes or low variation genes you should enable the checkbox in the `Data` section

**step 3)** Go to [interopt.ir](https://interopt.ir/). On the left panel in the `Data` section set Datasets as `Uploaded dataset`. Then Click on the Browse button and select your expression CT values file. leave other options on default mode and click on `Run Experiment`

**step 4)** After the process is finished, on the right panel you two tables are generated:

- the `Weighting Result` shows the weights and stability measures of each combination of the reference genes.
- the `Aggregated Genes` shows the aggregated reference genes after applying the weighted mean.

To download each of these tables in csv or excel format, click on the Download button at the left bottom corner of the tables. each of the rows in `Aggreagated Genes` table can be used like an artificial reference gene.

### UseCase 2: Selecting best combination of reference genes
The goal here is to find the best combination of reference genes from a set. For this case you should provide one of the following data with samples from your target biological condition.  

1. CT values from a qPCR array of large number of genes
2. CT values of 4 or more reference genes

The steps are exactly like previous section. After the processing is finished you can sort the `Weighting Result` table by clicking on the SD column to see the most stable combinations (with the lowest SD)


