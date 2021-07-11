import React, { useState } from 'react';
import DatePicker from 'react-datepicker';
import 'react-datepicker/dist/react-datepicker.css';
import './DatePickerInput.scss';
import datePickerIcon from '../../../../assets/images/date-picker-icon.svg';

const DatePickerInput = (props) => {
  const { label, value, onChange = () => {} } = props;

  const [date, setdate] = useState(value || new Date());

  const handleDateSelect = (date: Date) => {
    setdate(date);
    onChange(date);
  };

  const handleDateChange = (date: Date) => {
    setdate(date);
    onChange(date);
  };

  return (
    <div className="form-element">
      <label className="form__label">{label}</label>
      <div className="custom-container">
        <span></span>
        <DatePicker
          selected={date}
          onSelect={handleDateSelect} //when day is clicked
          onChange={handleDateChange} //only when value has changed
          showMonthDropdown={true}
          showYearDropdown={true}
        />
      </div>
      <div className="datepicker-icon">
        <img src={datePickerIcon}></img>
      </div>
    </div>
  );
};

export default DatePickerInput;
