import datetime
import math
import logging

def coerceDashedDate(ds):
    dt = datetime.datetime.strptime(ds, "%Y-%m-%d-%H-%M-%S")
    #year, month, day, hour, minute, second = [int(j) for j in ds.split("-")]
    #dt = datetime.datetime(year, month, day, hour=hour, minute=minute, second=second)
    return dt.timestamp()

def coerceDateTime(ds):
    #dt = datetime.datetime.strptime(ds, "%m/%d/%Y %H:%M:%S %p")
    # print(ds)
    date, time, noon = ds.split()
    month, day, year = [int(each) for each in date.split("/")]
    hour, minute, second = [int(each) for each in time.split(":")]
    if noon == "PM" and hour < 12:
        hour += 12
    elif noon == "AM" and hour == 12:
        hour = 0
    dt = datetime.datetime(year, month, day, hour=hour, minute=minute, second=second)
    return dt.timestamp()

def coerceStringDate(ds):
    dt = datetime.datetime.strptime(ds, "%d %b %Y %H:%M")
    return dt.timestamp()

def now(asstr = False, tz = datetime.timezone.utc):
    """
    Wrapper around datetime.now()

    A convenience function for returning the current time as a ISO 8601 or as a
    unix timestamp.
    """
    dt = datetime.datetime.now(tz = tz)
    if asstr:
        return dt.strftime("%Y-%m-%d %H:%M:%S")
    else:
        return dt.timestamp()

def _infer_timestamp_from(headers, spec = None, tz = datetime.timezone.utc):
    """
    Convenience function for timestamping

    Given a set of headers, and an optional specification, return an array
    containing column indices from which a timestamp in a given row can be
    computed, as well as the function which will compute the timestamp given
    the returned array.

    Parameters
    ----------
    headers : array
        An array of strings. If `spec` is not supplied, must contain either "uts"
        (float) or "timestep" (ISO 8601).

    spec : dict, optional
        A specification of timestamp elements with associated column indices and
        optional formats. Currently accepted combinations of keys are: "uts";
        "timestamp"; "date" and / or "time".
    
    tz : datetime.timezone, optional
        Timezone to use for conversion. By default, UTC is used.
    """
    if spec is not None:
        if "uts" in spec:
            return [spec["uts"]], float
        if "timestamp" in spec:
            if isinstance(spec["timestamp"], int):
                logging.debug("dateutils: Assuming specified column containing "
                              "the timestamp is in ISO 8601 format")
                def retfunc(value):
                    dtn = datetime.datetime.fromisoformat(value)
                    dt = datetime.datetime(year = dtn.year, month= dtn.month,
                                           day = dtn.day, hour = dtn.hour,
                                           minute = dtn.minute, second = dtn.second,
                                           tzinfo = tz)
                    return dt.timestamp()
                return [spec["timestamp"]], retfunc
            else:
                def retfunc(value):
                    dtn = datetime.datetime.strptime(value, spec["timestamp"][1])
                    dt = datetime.datetime(year = dtn.year, month= dtn.month,
                                           day = dtn.day, hour = dtn.hour,
                                           minute = dtn.minute, second = dtn.second,
                                           tzinfo = tz)
                    return dt.timestamp()
                return [spec["timestamp"][0]], retfunc
        if "date" in spec or "time" in spec:
            specdict = {
                "date": datetime.datetime.fromtimestamp(0, tz = tz).timestamp,
                "time": datetime.datetime.fromtimestamp(0, tz = tz).timestamp
            }
            cols = [None, None]
            if "date" in spec:
                if isinstance(spec["date"], int):
                    logging.debug("dateutils: Assuming specified column containing "
                                "the date is in ISO 8601 format")
                    def datefn(value):
                        dtn = datetime.datetime.fromisoformat(value)
                        dt = datetime.datetime(year = dtn.year, month= dtn.month,
                                           day = dtn.day, hour = dtn.hour,
                                           minute = dtn.minute, second = dtn.second,
                                           tzinfo = tz)
                        return dt.timestamp()
                    cols[0] = spec["date"]
                else:
                    def datefn(value):
                        dtn = datetime.datetime.strptime(value, spec["date"][1])
                        dt = datetime.datetime(year = dtn.year, month= dtn.month,
                                           day = dtn.day, hour = dtn.hour,
                                           minute = dtn.minute, second = dtn.second,
                                           tzinfo = tz)
                        return dt.timestamp()
                    cols[0] = spec["date"][0]
                specdict["date"] = datefn
            if "time" in spec:
                if isinstance(spec["time"], int):
                    logging.debug("dateutils: Assuming specified column containing "
                                  "the time is in ISO 8601 format")
                    def timefn(value):
                        t = datetime.time.fromisoformat(value)
                        td = datetime.timedelta(hours = t.hour, minutes = t.minute,
                                                seconds = t.second)
                        return td.total_seconds()
                    cols[1] = spec["time"]
                else:
                    def timefn(value):
                        t = datetime.datetime.strptime(value, spec["time"][1])
                        td = datetime.timedelta(hours = t.hour, minutes = t.minute, seconds = t.second)
                        return td.total_seconds()
                    cols[1] = spec["time"][0]
                specdict["time"] = timefn
            if cols[0] is None:
                return [cols[1]], specdict["time"]
            elif cols[1] is None:
                return [cols[0]], specdict["date"]
            else:
                def retfn(date, time):
                    return specdict["date"](date) + specdict["time"](time)
                return cols, retfn
    elif "uts" in headers:
        logging.info("dateutils: No timestamp spec provided, assuming column 'uts' "
                     "is a valid unix timestamp")
        return [headers.index("uts")], float
    elif "timestamp" in headers:
        logging.info("dateutils: No timestamp spec provided, assuming column 'timestamp' "
                     "is a valid ISO 8601 timestamp")
        def retfunc(value):
            dtn = datetime.datetime.fromisoformat(value)
            dt = datetime.datetime(year = dtn.year, month= dtn.month,
                                   day = dtn.day, hour = dtn.hour,
                                   minute = dtn.minute, second = dtn.second,
                                   tzinfo = tz)
            return dt.timestamp()
        return [headers.index("timestamp")], retfunc
    else:
        assert True, \
            logging.error("dateutils: A valid timestamp could not be deduced.")
