import React from 'react';
import './AddOption.scss';
import addOptionImg from '../../../../assets/images/add-option.svg';

const AddOption = (props) => {
  const { label = '', onClick = () => {} } = props;

  const handleClick = (event) => {
    event && event.stopPropagation && event.stopPropagation();
    onClick();
  };
  return (
    <div className="add-option" onClick={handleClick}>
      <img alt="Add" src={addOptionImg}></img>
      <a className="form__label-small">{label}</a>
    </div>
  );
};

export default AddOption;
