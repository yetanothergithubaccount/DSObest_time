# DSObest_time
Find out about the best time during a night to observe a list of deep sky objects.

This tool helps to produce a report in pdf-format of all the observable deep sky 
objects from the desired catalogue.

## Usage
```
python3 DSO_observation_planning.py --tonight --catalogue Messier --moon --configuration Frankfurt # check Messier DSO's for tonight, consider moon

python3 DSO_observation_planning.py --tonight --catalogue Caldwell --moon --configuration Frankfurt # check Caldwell DSO's for tonight, consider moon

python3 DSO_observation_planning.py --tonight --catalogue South --moon --configuration Windhoek # check some southern hemisphere  DSO's for tonight, consider moon
```
## Result
The resulting PDF-document for a list of well-observable deep sky objects above Frankfurt produced with
```
python3 DSO_observation_planning.py --tonight --catalogue Messier --moon --configuration Frankfurt # check Messier DSO's for tonight, consider moon
```
may look like this:

![Screenshot of PDF-report](https://github.com/yetanothergithubaccount/DSObest_time/blob/main/pdfreportscreenshot.png)
