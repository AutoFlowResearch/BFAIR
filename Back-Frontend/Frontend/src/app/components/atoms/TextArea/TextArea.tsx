import React, { ChangeEvent, useState } from 'react';
import './TextArea.style.scss';

const TextArea = (props) => {
  const { label = '', value = '', onChange = () => {} } = props;
  const [inputValue, setInputValue] = useState(value);
  const handleChange = (event: ChangeEvent<HTMLTextAreaElement>) => {
    const newVal = event.target?.value;
    setInputValue(newVal);
    onChange(newVal);
  };
  return (
    <div className='form-element textarea-container'>
      <label className='form__label'>{label}</label>
      <textarea
        className='form-textarea'
        onChange={handleChange}
        placeholder={`Enter a ${label}...`}
        value={inputValue}
      ></textarea>
    </div>
  );
};

export default TextArea;
